// Mitogenome extraction with GetOrganelle.
//
// GetOrganelle ships a per-organelle "label" database (seed sequences + HMM
// profiles). It is downloaded on first use into ~/.GetOrganelle, which is
// ephemeral inside a container - so we materialise the DB as a Nextflow
// channel: MITO_DB downloads it once into a directory that downstream
// MITOGENOME calls mount as input. Set params.organelle_db to skip the
// download (point at an existing GetOrganelle library directory).
//
// Seeding: GetOrganelle recruits its starting reads by bowtie2 against a seed
// FASTA (default: the bundled seed DB for organelle_type). For divergent taxa
// the bundled animal_mt seeds can recruit no true mitochondrial reads at all;
// the coverage estimate is then taken from nuclear hits and the run fails.
// meta.mito_seed replaces the default seed with a close relative's sequence
// (same genus or closer; a COX1 barcode is enough). It comes from the
// samplesheet column, from --mito_seed_library (then meta.mito_seed_records
// names the records of that sample, chosen in main.nf by species, sample id or
// genus), or from --mito_seed. Gap characters are removed, so aligned FASTAs
// work as seeds. The label DB used to judge the final graph stays the default.

include { FIND_SEED } from './mito_seed.nf'

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
    // Keeping all reads (params.mito_keep_all_reads) is slower on deep
    // libraries, so the walltime also scales with the retry attempt.
    time { 4.h * task.attempt }
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
    tuple val(meta),  path(reads), path(seed, stageAs: 'mito_seed/*'), val(seed_src)
    tuple val(otype), path(db)

    output:
    tuple val(meta), path("${meta.id}.mito.fasta"),     emit: fasta, optional: true
    tuple val(meta), path("${meta.id}.mito.log.txt"),   emit: logfile
    tuple val(meta), path("${meta.id}.getorganelle_mqc.tsv"), emit: summary

    script:
    def in_args   = meta.single_end ? "-u ${reads[0]}" : "-1 ${reads[0]} -2 ${reads[1]}"
    // seed_src: 'custom' / 'library' (meta.mito_seed), 'self' (FIND_SEED
    // contigs from the sample's own assembly) or 'default' (no seed file).
    def use_seed  = seed_src != 'default'
    def recs      = seed_src == 'library' ? (meta.mito_seed_records ?: []) : []
    def seed_arg  = use_seed ? '-s seed_clean.fasta' : ''
    def seed_name = seed_src == 'default' ? 'default'
                  : seed_src == 'self' ? 'self:own_assembly'
                  : recs ? "library:${recs.size() <= 2 ? recs.join('+') : recs[0] + '+' + (recs.size() - 1) + '_more'}"
                  : "custom:${file(meta.mito_seed).name}"
    // Keep the selected records (all when none are named) and drop alignment
    // gap characters, which bowtie2 cannot index.
    def seed_prep = use_seed ? """
    awk -v want='${recs.join(' ')}' '
        BEGIN { n = split(want, w, " "); for (i = 1; i <= n; i++) keep[w[i]] = 1 }
        /^>/  { split(substr(\$0, 2), h, /[ \t]+/); p = (n == 0 || (h[1] in keep)); if (p) print; next }
        p     { gsub(/[-.[:space:]]/, ""); if (length(\$0)) print }
    ' ${seed} > seed_clean.fasta
    if ! grep -q '^[^>]' seed_clean.fasta; then
        echo "mito seed: no sequence selected from ${seed} (records: ${recs.join(' ') ?: 'all'})" >&2
        exit 1
    fi
""" : ''
    // GetOrganelle subsamples reads to a target organelle coverage estimated
    // from the seed reads. A bad seed inflates that estimate (950x to 2471x on
    // the pilargid test set, true value ~40x) and most mito reads are dropped,
    // so by default every read is kept.
    def keep_arg  = params.mito_keep_all_reads ? '--reduce-reads-for-coverage inf --max-reads inf' : ''
    """
    set -euo pipefail
    ${seed_prep}
    set +e
    GETORG_PATH=\$PWD/${db} \\
    get_organelle_from_reads.py \\
        ${in_args} \\
        -F ${otype} \\
        ${seed_arg} \\
        ${keep_arg} \\
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
Sample	Status	Sequences	Total_bp	Longest_bp	Topology	Path_FASTAs	Seed
${meta.id}	\$status_label	\$seq_count	\$total_bp	\$longest_bp	\$topology	\$path_count	${seed_name}
EOF

    exit 0
    """
}

workflow MITO {
    take:
    reads_ch         // channel: tuple(meta, reads)
    contigs_ch       // channel: tuple(meta, contigs), the MEGAHIT assemblies (self-seeding)

    main:
    if (params.organelle_db) {
        // User supplied a pre-downloaded DB.
        db_ch = Channel.value(tuple(params.organelle_type, file(params.organelle_db, checkIfExists: true)))
    } else {
        db_req_ch = reads_ch.map { meta, reads -> params.organelle_type }.first()
        MITO_DB(db_req_ch)
        db_ch = MITO_DB.out
    }

    // Seed per sample. A given seed (samplesheet, library, --mito_seed) is
    // used directly; otherwise, with params.mito_self_seed, the sample waits
    // for its assembly and FIND_SEED looks for its own mito contigs; else the
    // default seeds ([] stages nothing).
    split_ch = reads_ch.branch { meta, reads ->
        given: meta.mito_seed
        self:  params.mito_self_seed
        plain: true
    }
    given_ch = split_ch.given.map { meta, reads ->
        tuple(meta, reads, file(meta.mito_seed), meta.mito_seed_records ? 'library' : 'custom')
    }
    plain_ch = split_ch.plain.map { meta, reads -> tuple(meta, reads, [], 'default') }

    if (params.mito_self_seed) {
        ref_ch = Channel.value(file(params.mito_seed_ref, checkIfExists: true))
        FIND_SEED(contigs_ch.join(split_ch.self).map { meta, contigs, reads -> tuple(meta, contigs) }, ref_ch)
        self_ch = split_ch.self
            .join(FIND_SEED.out.seed, remainder: true)
            .filter { meta, reads, seed -> reads != null }
            .map { meta, reads, seed -> seed ? tuple(meta, reads, seed, 'self') : tuple(meta, reads, [], 'default') }
        seed_summary = FIND_SEED.out.summary
    } else {
        self_ch      = Channel.empty()
        seed_summary = Channel.empty()
    }

    MITOGENOME(given_ch.mix(self_ch, plain_ch), db_ch)

    emit:
    fasta        = MITOGENOME.out.fasta
    logfile      = MITOGENOME.out.logfile
    summary      = MITOGENOME.out.summary
    seed_summary = seed_summary
}
