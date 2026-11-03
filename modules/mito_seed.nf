// Self-seeding for GetOrganelle: find the sample's own mitochondrial contigs
// in its MEGAHIT assembly and use them as the seed.
//
// Approach after the seed extraction Mateo developed for the syllid skims
// (profile search of each sample's own assembly, own mito recognised by its
// read depth). Here: tblastn with the 13 mitochondrial protein-coding genes of
// eight annelids (assets/mito_pcg_annelida.faa; override with
// --mito_seed_ref) under the mitochondrial genetic code; bin/find_mito_seed.py
// drops hits with stop codons (NUMTs) or little query coverage and keeps the
// deepest clean contigs. Only for samples without a seed from the samplesheet,
// --mito_seed_library or --mito_seed. No hit => GetOrganelle's default seeds.
//
// Runs in the GetOrganelle image, which ships BLAST+ and Python.

process FIND_SEED {
    tag "${meta.id}"
    container 'quay.io/biocontainers/getorganelle:1.7.7.1--pyhdfd78af_0'
    publishDir "${params.outdir}/mitogenome", mode: 'copy', saveAs: { fn -> "${meta.id}/${fn}" }

    cpus 4
    memory '8 GB'
    time '4h'
    // A failed search only means "no self seed": the sample then falls back
    // to the default seeds in MITO.
    errorStrategy 'ignore'

    input:
    tuple val(meta), path(contigs)
    path ref

    output:
    tuple val(meta), path("${meta.id}.self_seed.fasta"), emit: seed, optional: true
    tuple val(meta), path("${meta.id}.self_seed_mqc.tsv"), emit: summary

    script:
    """
    set -euo pipefail

    # Short contigs cannot hold a useful seed and dominate the search time.
    awk -v min=${params.mito_seed_min_contig} '
        /^>/ { if (h != "" && length(s) >= min) print h "\\n" s; h = \$0; s = ""; next }
        { s = s \$0 }
        END  { if (h != "" && length(s) >= min) print h "\\n" s }
    ' ${contigs} > long_contigs.fa

    if [ -s long_contigs.fa ]; then
        makeblastdb -in long_contigs.fa -dbtype nucl -out ctg > /dev/null
        tblastn -query ${ref} -db ctg \\
            -db_gencode ${params.mito_seed_gencode} \\
            -evalue 1e-10 -max_target_seqs 50 \\
            -outfmt '6 qseqid sseqid pident length qstart qend sstart send evalue bitscore qlen sseq' \\
            -num_threads ${task.cpus} > hits.tsv
    else
        : > hits.tsv
    fi

    find_mito_seed.py \\
        --hits hits.tsv \\
        --contigs long_contigs.fa \\
        --sample ${meta.id} \\
        --depth-frac ${params.mito_seed_depth_frac} \\
        --out-fasta ${meta.id}.self_seed.fasta \\
        --out-tsv ${meta.id}.self_seed_mqc.tsv
    """
}
