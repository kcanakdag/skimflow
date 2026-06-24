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

    cpus 8
    // GetOrganelle runs SPAdes internally; on full-size libraries a fixed 6 GB
    // gets OOM-killed. Start at 8 cores x 2 GB/core = 16 GB (the scc-cpu node
    // ratio) and scale with attempt (16 -> 32 -> 48 -> 64 GB) so the GWDG
    // OOM-retry recovers instead of dropping the mitogenome.
    memory { 16.GB * task.attempt }
    time '4h'
    maxRetries 3

    // GetOrganelle can return non-zero when it cannot finish the assembly
    // (e.g. too few mitochondrial reads). Capture that in a summary instead
    // of taking down the rest of the run.
    //
    // On OOM (exit 137, see the kill-signal guard in the script) retry with
    // the next memory tier; once retries are spent, ignore so one stubborn
    // sample can't abort the whole batch.
    errorStrategy { task.attempt <= 3 && (task.exitStatus in 130..143) ? 'retry' : 'ignore' }

    input:
    tuple val(meta),  path(reads)
    tuple val(otype), path(db)

    output:
    tuple val(meta), path("${meta.id}.mito.fasta"),     emit: fasta, optional: true
    tuple val(meta), path("${meta.id}.mito.log.txt"),   emit: logfile
    tuple val(meta), path("${meta.id}.getorganelle_mqc.tsv"), emit: summary

    script:
    def in_args = meta.single_end ? "-u ${reads[0]}" : "-1 ${reads[0]} -2 ${reads[1]}"
    """
    set -euo pipefail

    set +e
    GETORG_PATH=\$PWD/${db} \\
    get_organelle_from_reads.py \\
        ${in_args} \\
        -F ${otype} \\
        -o mito_out \\
        -k 21,45,65,85,105 \\
        -t ${task.cpus} \\
        -R 10
    status=\$?
    set -e

    # If GetOrganelle was killed by a signal (e.g. the OOM-killer: SIGKILL ->
    # exit 137), propagate it so Nextflow's memory-scaling retry kicks in.
    # Without this the 'exit 0' below would mask the OOM as an empty success,
    # which also poisons -resume (a no_fasta result gets cached as done).
    # Genuine GetOrganelle errors (too few reads -> small exit codes) fall
    # through and are tolerated as an empty result.
    if [ "\$status" -ge 128 ]; then
        echo "GetOrganelle killed by signal (exit \$status), likely OOM; failing task to trigger retry" >&2
        exit "\$status"
    fi

    if [ -f mito_out/get_org.log.txt ]; then
        cp mito_out/get_org.log.txt ${meta.id}.mito.log.txt
    else
        printf "GetOrganelle exited with status %s before writing get_org.log.txt\\n" "\$status" > ${meta.id}.mito.log.txt
    fi

    # GetOrganelle's documented final FASTA output is path_sequence.fasta.
    # Multiple files can represent alternative graph paths; annotate the first
    # deterministically and record the count in the report summary.
    if [ -d mito_out ]; then
        path_count=\$(find mito_out -maxdepth 1 -type f -name '*path_sequence*.fasta' | sort | wc -l | awk '{print \$1}')
        found=\$(find mito_out -maxdepth 1 -type f -name '*path_sequence*.fasta' | sort | head -1 || true)
    else
        path_count=0
        found=""
    fi
    if [ -n "\$found" ]; then
        cp "\$found" ${meta.id}.mito.fasta
    fi

    if [ -f ${meta.id}.mito.fasta ]; then
        status_label="ok"
        seq_count=\$(grep -c '^>' ${meta.id}.mito.fasta || true)
        total_bp=\$(awk 'BEGIN{n=0} /^>/{next} {gsub(/[[:space:]]/, ""); n+=length(\$0)} END{print n+0}' ${meta.id}.mito.fasta)
        longest_bp=\$(awk 'BEGIN{m=0;n=0} /^>/{if(n>m)m=n; n=0; next} {gsub(/[[:space:]]/, ""); n+=length(\$0)} END{if(n>m)m=n; print m+0}' ${meta.id}.mito.fasta)
        topology=\$(grep '^>' ${meta.id}.mito.fasta | grep -qi circular && printf circular || printf unknown)
    else
        status_label="no_fasta"
        seq_count=0
        total_bp=0
        longest_bp=0
        topology="NA"
    fi

cat > ${meta.id}.getorganelle_mqc.tsv <<EOF
# id: "getorganelle_mitogenome"
# section_name: "Mitogenome assembly"
# description: "GetOrganelle mitochondrial FASTA summary."
# plot_type: "table"
Sample	Status	Sequences	Total_bp	Longest_bp	Topology	Path_FASTAs
${meta.id}	\$status_label	\$seq_count	\$total_bp	\$longest_bp	\$topology	\$path_count
EOF

    exit 0
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
        db_req_ch = reads_ch.map { meta, reads -> params.organelle_type }.first()
        MITO_DB(db_req_ch)
        db_ch = MITO_DB.out.collect(flat: false).map { it[0] }
    }

    MITOGENOME(reads_ch, db_ch)

    emit:
    fasta   = MITOGENOME.out.fasta
    logfile = MITOGENOME.out.logfile
    summary = MITOGENOME.out.summary
}
