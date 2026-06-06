// Long-read QC: Filtlong length/quality filtering. Runs only for long-read
// samples. Filtlong has no output flag; it streams filtered reads to stdout,
// so we pipe to gzip.

process LR_QC {
    tag "${meta.id}"
    container 'quay.io/biocontainers/filtlong:0.3.1--h077b44d_0'
    publishDir "${params.outdir}/long_read_qc", mode: 'copy', saveAs: { fn -> "${meta.id}/${fn}" }

    cpus 1
    memory '4 GB'
    time '1h'

    input:
    tuple val(meta), path(long_reads)

    output:
    tuple val(meta), path("${meta.id}.filtlong.fastq.gz"), emit: reads

    script:
    """
    filtlong \\
        --min_length ${params.filtlong_min_length} \\
        --keep_percent ${params.filtlong_keep_percent} \\
        ${long_reads} | gzip > ${meta.id}.filtlong.fastq.gz
    """
}
