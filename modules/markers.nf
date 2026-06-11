// Marker extraction with BUSCO.
// Runs in genome mode against the assembly (MEGAHIT for short reads, Flye for
// long reads); the BUSCO short summary gives the standard C/F/D/M numbers used
// in skim assembly papers.

process BUSCO_DB {
    tag "${lineage}"
    container 'quay.io/biocontainers/busco:6.0.0--pyhdfd78af_3'
    publishDir "${params.outdir}/markers", mode: 'copy'

    cpus 1
    memory '2 GB'

    input:
    val lineage

    output:
    tuple val(lineage), path('busco_downloads')

    script:
    """
    mkdir -p busco_downloads
    busco --download_path busco_downloads --download ${lineage}
    """
}

process MARKERS {
    tag "${meta.id}"
    container 'quay.io/biocontainers/busco:6.0.0--pyhdfd78af_3'
    publishDir "${params.outdir}/markers", mode: 'copy', saveAs: { fn -> "${meta.id}/${fn}" }

    cpus 8
    memory '12 GB'
    time '4h'

    // BUSCO can fail when a tiny test assembly has zero hits. Don't take
    // the rest of the pipeline down with it.
    errorStrategy 'ignore'

    input:
    tuple val(meta),    path(contigs)
    tuple val(lineage), path(busco_dl)

    output:
    tuple val(meta), path("${meta.id}_busco/short_summary*.txt"), emit: summary, optional: true
    tuple val(meta), path("${meta.id}_busco"),                    emit: dir,     optional: true

    script:
    """
    busco \\
        --in ${contigs} \\
        --mode genome \\
        --lineage_dataset ${lineage} \\
        --download_path ${busco_dl} \\
        --offline \\
        --cpu ${task.cpus} \\
        --out ${meta.id}_busco \\
        --out_path .
    """
}

workflow MARKER_EXTRACTION {
    take:
    contigs_ch       // channel: tuple(meta, contigs)

    main:
    if (params.busco_db) {
        db_ch = Channel.value(tuple(params.busco_lineage, file(params.busco_db, checkIfExists: true)))
    } else {
        db_req_ch = contigs_ch.map { meta, contigs -> params.busco_lineage }.first()
        BUSCO_DB(db_req_ch)
        db_ch = BUSCO_DB.out.collect(flat: false).map { it[0] }
    }

    MARKERS(contigs_ch, db_ch)

    emit:
    summary = MARKERS.out.summary
    dir     = MARKERS.out.dir
}
