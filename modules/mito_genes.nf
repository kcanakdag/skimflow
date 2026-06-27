// Harvest per-gene sequences from each sample's MITOS2 annotation and
// aggregate them across all samples into one multi-FASTA per gene, plus a
// gene-occupancy matrix. Parsing logic lives in bin/harvest_mito_genes.py and
// bin/aggregate_mito_genes.py (auto-added to PATH by Nextflow).

process GENE_HARVEST {
    tag "${meta.id}"
    container params.gene_container
    publishDir "${params.outdir}/genes/per_sample", mode: 'copy', saveAs: { fn -> "${meta.id}/${fn}" }

    cpus 1
    memory '2 GB'
    time '30m'
    // Same policy as the annotation steps: a sample MITOS2 could not annotate
    // simply contributes nothing rather than aborting the batch.
    errorStrategy 'ignore'

    input:
    tuple val(meta), path(mitos2_dir)

    output:
    path 'seqs/*',                            emit: genes,     optional: true
    path "${meta.id}.occupancy.tsv",          emit: occupancy, optional: true

    script:
    def species = meta.species ?: meta.id
    """
    harvest_mito_genes.py \\
        --mitos-dir ${mitos2_dir} \\
        --sample-id ${meta.id} \\
        --species-id '${species}' \\
        --outdir .
    """
}

process GENE_AGGREGATE {
    tag 'mito_genes'
    container params.gene_container
    publishDir "${params.outdir}/genes", mode: 'copy'

    cpus 1
    memory '2 GB'
    time '30m'

    input:
    path seqs, stageAs: 'seqs/*'
    path occ,  stageAs: 'occ/*'

    output:
    path 'nt/*.fasta',                     emit: nt,      optional: true
    path 'aa/*.faa',                       emit: aa,      optional: true
    path 'occupancy.tsv',                  emit: matrix
    path 'mito_genes_occupancy_mqc.tsv',   emit: summary

    script:
    """
    mkdir -p seqs occ
    aggregate_mito_genes.py --seqs-dir seqs --occ-dir occ --outdir .
    """
}

workflow MITO_GENES {
    take:
    mitos2_dir_ch    // channel: tuple(meta, mitos2_dir)

    main:
    GENE_HARVEST(mitos2_dir_ch)

    seqs_ch = GENE_HARVEST.out.genes.collect().ifEmpty([])
    occ_ch  = GENE_HARVEST.out.occupancy.collect().ifEmpty([])
    GENE_AGGREGATE(seqs_ch, occ_ch)

    emit:
    summary = GENE_AGGREGATE.out.summary
    nt      = GENE_AGGREGATE.out.nt
    aa      = GENE_AGGREGATE.out.aa
}
