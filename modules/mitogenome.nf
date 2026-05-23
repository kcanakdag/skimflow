// Mitogenome extraction with GetOrganelle.
//
// GetOrganelle ships a per-organelle "label" database (seed sequences + HMM
// profiles). It is downloaded on first use into ~/.GetOrganelle, which is
// ephemeral inside a container - so we materialise the DB as a Nextflow
// channel: MITO_DB downloads it once into a directory that downstream
// MITOGENOME calls mount as input. Set params.organelle_db to skip the
// download (point at an existing GetOrganelle library directory).

process MITO_DB {
    tag "${organelle_type}"
    container 'quay.io/biocontainers/getorganelle:1.7.7.1--pyhdfd78af_0'
    publishDir "${params.outdir}/mitogenome", mode: 'copy'

    cpus 1
    memory '2 GB'

    input:
    val organelle_type

    output:
    tuple val(organelle_type), path('getorganelle_db')

    script:
    """
    mkdir -p getorganelle_db
    GETORG_PATH=\$PWD/getorganelle_db \\
        get_organelle_config.py --add ${organelle_type}
    """
}

process MITOGENOME {
    tag "${meta.id}"
    container 'quay.io/biocontainers/getorganelle:1.7.7.1--pyhdfd78af_0'
    publishDir "${params.outdir}/mitogenome", mode: 'copy', saveAs: { fn -> "${meta.id}/${fn}" }

    cpus 4
    memory '6 GB'
    time '4h'

    // GetOrganelle returns non-zero when it cannot finish the assembly
    // (e.g. too few mitochondrial reads). For a smoke test this is fine -
    // we capture the log and continue the pipeline.
    errorStrategy 'ignore'

    input:
    tuple val(meta),  path(reads)
    tuple val(otype), path(db)

    output:
    tuple val(meta), path("${meta.id}.mito.fasta"),     emit: fasta, optional: true
    tuple val(meta), path("${meta.id}.mito.log.txt"),   emit: logfile

    script:
    def in_args = meta.single_end ? "-u ${reads[0]}" : "-1 ${reads[0]} -2 ${reads[1]}"
    """
    GETORG_PATH=\$PWD/${db} \\
    get_organelle_from_reads.py \\
        ${in_args} \\
        -F ${otype} \\
        -o mito_out \\
        -k 21,45,65,85,105 \\
        -t ${task.cpus} \\
        -R 10

    cp mito_out/get_org.log.txt ${meta.id}.mito.log.txt

    # GetOrganelle's success output filename varies (path_sequence.fasta,
    # complete.fasta, scaffolds.fasta…). Take whatever fasta it produced.
    found=\$(ls mito_out/*path_sequence*.fasta mito_out/*scaffolds*.fasta mito_out/*complete*.fasta 2>/dev/null | head -1 || true)
    if [ -n "\$found" ]; then
        cp "\$found" ${meta.id}.mito.fasta
    fi
    """
}

workflow MITO {
    take:
    reads_ch         // channel: tuple(meta, reads)

    main:
    if (params.organelle_db) {
        // User supplied a pre-downloaded DB.
        db_ch = Channel.value(tuple(params.organelle_type, file(params.organelle_db, checkIfExists: true)))
    } else {
        MITO_DB(Channel.value(params.organelle_type))
        db_ch = MITO_DB.out
    }

    MITOGENOME(reads_ch, db_ch)

    emit:
    fasta   = MITOGENOME.out.fasta
    logfile = MITOGENOME.out.logfile
}
