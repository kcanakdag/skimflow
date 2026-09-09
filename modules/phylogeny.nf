// Phylogenetics on the harvested per-gene FASTAs from MITO_GENES.
//
// Per-gene alignment then concatenation is the standard treatment for
// mitochondrial genes: whole-mitogenome alignment breaks on gene rearrangements
// (notable in annelids), the unalignable control region, and the arbitrary
// linearisation point of a circular molecule. So each gene is aligned (MAFFT)
// and trimmed (trimAl) on its own, an optional per-gene tree is built, and the
// trimmed alignments are concatenated into a partitioned supermatrix for one
// species tree (IQ-TREE). Both sequence types run through the same path:
//   nt -> 13 PCG + 2 rRNA (DNA)     aa -> 13 PCG (protein)
//
// Fail-soft like the rest of the pipeline: any single-gene or supermatrix tree
// that cannot be built (too few taxa, tool error) is skipped with a note rather
// than aborting the run. Concatenation logic lives in bin/concat_alignments.py.

process ALIGN {
    tag "${seqtype}:${gene}"
    container params.mafft_container
    publishDir "${params.outdir}/phylogeny", mode: 'copy', saveAs: { fn -> "${seqtype}/alignments/${fn}" }

    cpus 2
    memory '2 GB'
    time '1h'
    errorStrategy 'ignore'

    input:
    tuple val(seqtype), val(gene), path(fasta)

    output:
    tuple val(seqtype), val(gene), path("${gene}.aln.fasta"), emit: aln

    script:
    """
    set -euo pipefail
    n=\$(grep -c '^>' ${fasta} || true)
    if [ "\$n" -lt 2 ]; then
        # MAFFT needs >=2 sequences; a single-taxon gene is trivially aligned.
        cp ${fasta} ${gene}.aln.fasta
    else
        mafft ${params.mafft_args} ${fasta} > ${gene}.aln.fasta
    fi
    """
}

process TRIM {
    tag "${seqtype}:${gene}"
    container params.trimal_container
    publishDir "${params.outdir}/phylogeny", mode: 'copy', saveAs: { fn -> "${seqtype}/trimmed/${fn}" }

    cpus 1
    memory '2 GB'
    time '30m'
    errorStrategy 'ignore'

    input:
    tuple val(seqtype), val(gene), path(aln)

    output:
    tuple val(seqtype), val(gene), path("${gene}.trimmed.fasta"), emit: trimmed

    script:
    """
    set -euo pipefail
    # trimAl can strip every column on very gappy/divergent input; if it errors
    # or yields an empty file, keep the untrimmed alignment so the gene still
    # contributes to the supermatrix.
    if trimal -in ${aln} -out ${gene}.trimmed.fasta -${params.trimal_mode} 2>trimal.err && [ -s ${gene}.trimmed.fasta ]; then
        :
    else
        cp ${aln} ${gene}.trimmed.fasta
    fi
    """
}

process GENE_TREE {
    tag "${seqtype}:${gene}"
    container params.iqtree_container
    publishDir "${params.outdir}/phylogeny", mode: 'copy', saveAs: { fn -> "${seqtype}/gene_trees/${fn}" }

    cpus 2
    memory '4 GB'
    time '2h'
    errorStrategy 'ignore'

    input:
    tuple val(seqtype), val(gene), path(aln)

    output:
    path "${gene}.treefile",    emit: tree, optional: true
    path "${gene}.iqtree",      optional: true
    path "${gene}.contree",     optional: true
    path "${gene}.log",         optional: true
    path "${gene}.skipped.txt", optional: true

    script:
    """
    set -euo pipefail
    n=\$(grep -c '^>' ${aln} || true)
    if [ "\$n" -lt 4 ]; then
        # A single-gene ML tree with UFBoot needs at least 4 taxa.
        echo "gene ${gene} (${seqtype}): only \$n taxa (<4); single-gene tree skipped" > ${gene}.skipped.txt
        exit 0
    fi
    iqtree2 \\
        -s ${aln} \\
        -m MFP \\
        -B ${params.iqtree_bootstrap} \\
        -alrt ${params.iqtree_alrt} \\
        -T AUTO -ntmax ${task.cpus} \\
        --prefix ${gene} -redo
    """
}

process CONCAT {
    tag "${seqtype}"
    container params.gene_container
    publishDir "${params.outdir}/phylogeny", mode: 'copy', saveAs: { fn -> "${seqtype}/${fn}" }

    cpus 1
    memory '2 GB'
    time '30m'

    input:
    tuple val(seqtype), path(alns, stageAs: 'aln/*')

    output:
    tuple val(seqtype), path('supermatrix.fasta'), path('partitions.nex'), emit: supermat
    path 'concat_stats.tsv', emit: stats

    script:
    """
    concat_alignments.py --aln-dir aln --outdir . --min-taxa ${params.phylo_min_taxa_per_gene}
    """
}

process SUPER_TREE {
    tag "${seqtype}"
    container params.iqtree_container
    publishDir "${params.outdir}/phylogeny", mode: 'copy', saveAs: { fn -> "${seqtype}/${fn}" }

    cpus 4
    memory '8 GB'
    time '6h'
    errorStrategy 'ignore'

    input:
    tuple val(seqtype), path(supermatrix), path(partitions)

    output:
    path "${seqtype}_species.*",         emit: tree, optional: true
    tuple val(seqtype), path("${seqtype}_species.treefile"), emit: treefile, optional: true
    path "${seqtype}.phylo_mqc.tsv",     emit: summary

    script:
    """
    set -euo pipefail
    ntax=\$(grep -c '^>' ${supermatrix} 2>/dev/null || echo 0)
    ngenes=\$(grep -c 'charset' ${partitions} 2>/dev/null || echo 0)
    alnlen=\$(awk 'NR==2{print length(\$0); exit}' ${supermatrix} 2>/dev/null || echo 0)
    alnlen=\${alnlen:-0}

    if [ "\$ntax" -lt 4 ]; then
        status="skipped_lt4taxa"
    else
        set +e
        iqtree2 \\
            -s ${supermatrix} \\
            -p ${partitions} \\
            -m MFP \\
            -B ${params.iqtree_bootstrap} \\
            -alrt ${params.iqtree_alrt} \\
            -T AUTO -ntmax ${task.cpus} \\
            --prefix ${seqtype}_species -redo > ${seqtype}_species.run.log 2>&1
        set -e
        if [ -f ${seqtype}_species.treefile ]; then status="ok"; else status="failed"; fi
    fi

cat > ${seqtype}.phylo_mqc.tsv <<EOF
# id: "phylogeny_${seqtype}"
# section_name: "Phylogeny (${seqtype})"
# description: "Supermatrix species tree from MAFFT + trimAl + IQ-TREE over harvested mito genes (${seqtype}). Drawn below; the .treefile is also published for FigTree/iTOL."
# plot_type: "table"
Seq_type	Status	Taxa	Genes	Alignment_length
${seqtype}	\$status	\$ntax	\$ngenes	\$alnlen
EOF

    exit 0
    """
}

process TREE_PLOT {
    tag "${seqtype}"
    container params.gene_container
    publishDir "${params.outdir}/phylogeny", mode: 'copy', saveAs: { fn -> "${seqtype}/${fn}" }

    cpus 1
    memory '1 GB'
    time '15m'
    errorStrategy 'ignore'

    input:
    tuple val(seqtype), path(treefile)

    output:
    path "${seqtype}_tree_mqc.html", emit: html

    script:
    """
    tree_to_svg.py --treefile ${treefile} --seqtype ${seqtype} --out ${seqtype}_tree_mqc.html
    """
}

workflow PHYLOGENY {
    take:
    nt_ch    // channel: MITO_GENES.out.nt  (per-gene nucleotide FASTAs)
    aa_ch    // channel: MITO_GENES.out.aa  (per-gene amino-acid FASTAs)

    main:
    // Fan each per-gene FASTA out into (seqtype, gene_stem, file). baseName of
    // COX1.fasta / COX1.faa is 'COX1', which is the gene stem.
    nt_genes = nt_ch.flatten().map { f -> tuple('nt', f.baseName, f) }
    aa_genes = aa_ch.flatten().map { f -> tuple('aa', f.baseName, f) }
    genes    = nt_genes.mix(aa_genes)

    ALIGN(genes)
    TRIM(ALIGN.out.aln)

    if (params.phylo_gene_trees) {
        GENE_TREE(TRIM.out.trimmed)
    }

    // Group trimmed alignments by sequence type -> one supermatrix per type.
    concat_in = TRIM.out.trimmed.map { st, g, f -> tuple(st, f) }.groupTuple()
    CONCAT(concat_in)
    SUPER_TREE(CONCAT.out.supermat)

    // Render each species tree that was actually built into an inline-SVG
    // MultiQC section. SUPER_TREE emits treefile only when a tree exists, so
    // skipped/failed seqtypes simply produce no plot.
    TREE_PLOT(SUPER_TREE.out.treefile)

    emit:
    summary  = SUPER_TREE.out.summary
    treeplot = TREE_PLOT.out.html
}
