// Genome-size + coverage estimation with RESPECT.
//
// RESPECT (Sayyari et al. 2022, Bioinformatics) is a k-mer-based estimator
// designed specifically for low-coverage skim data. From a single FASTQ
// pair it produces:
//   - genome_length         : haploid genome size (bp)
//   - sequencing_coverage   : per-base coverage (×)
//   - HCRM repeat profile   : high-copy repeat fraction (used internally)
// The "is the data good enough?" decision is made on these two columns.
//
// RESPECT requires the Gurobi optimization solver (academic licence is
// free via Gurobi WLS but not on biocontainers), so we run it inside a
// custom container. See containers/respect/Containerfile for the build
// recipe. The Gurobi licence file is supplied at runtime via:
//   - params.gurobi_lic    -> bound into the container at /opt/gurobi/gurobi.lic
//   - or env GRB_LICENSE_FILE pointing at the same path on the host.
//
// If params.expected_genome_size_bp is set in the samplesheet, RESPECT
// is still run (we want its coverage estimate), but the override is
// surfaced alongside in the MultiQC report as a sanity check.

process RESPECT {
    tag "${meta.id}"
    container "${params.respect_container ?: 'local/respect:0.2'}"
    publishDir "${params.outdir}/genome_size", mode: 'copy', saveAs: { fn -> "${meta.id}/${fn}" }

    cpus 8
    // Full-size skim libraries (multi-GB) overflow a small heap during k-mer
    // counting and get OOM-killed. Start at 32 GB (8 cores x 4 GB/core on the
    // scc-cpu node) and scale with attempt (32 -> 64 -> 96 -> 128 GB). The
    // script below re-exits 137 on a swallowed count OOM so this scaling
    // actually engages via the GWDG retry (exit 130..143); without that guard
    // RESPECT's plain exit-1 masks the OOM and the retry never fires.
    memory { 32.GB * task.attempt }
    time '4h'
    maxRetries 3

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("respect_out/estimated-parameters.txt"), emit: estimates
    tuple val(meta), path("${meta.id}.genome_size.txt"),           emit: size
    tuple val(meta), path("${meta.id}.coverage.txt"),              emit: coverage
    tuple val(meta), path("${meta.id}_respect_mqc.tsv"),           emit: summary
    tuple val(meta), path("respect_out"),                          emit: dir

    script:
    def reads_args = meta.single_end ? "${reads[0]}" : "${reads[0]} ${reads[1]}"
    def expected   = meta.expected_size ? "${meta.expected_size}" : 'NA'
    """
    mkdir -p respect_in
    # Stage the reads under a single directory so we can hand RESPECT a -d
    # argument (it scans for FASTQ inputs there).
    for f in ${reads_args}; do
        ln -s "\$(readlink -f \$f)" respect_in/\$(basename \$f)
    done

    # RESPECT CLI:
    #  -d <paths>      scan dirs for FASTQ inputs
    #  -o <path>       output dir
    #  --decomp gzip   use Python gzip (explicit output path). The default
    #                  cl-gzip follows symlinks and writes the decompressed
    #                  file next to the symlink TARGET (in some other work
    #                  dir), then RESPECT can't clean it up. Python decomp
    #                  always writes alongside the .gz, which is what we
    #                  want with Nextflow staging.
    #  --threads N     RESPECT-level worker count (one per input file).
    # Run RESPECT capturing its full output. RESPECT swallows a fatal jellyfish
    # k-mer-count failure (e.g. an OOM-kill of the count subprocess): it logs
    # "it's skipped", drops the input, then exits 1 once every input was skipped
    # ("Number of processes must be at least 1"). That plain exit-1 hides the
    # underlying OOM, so the GWDG memory-scaling retry (which only fires on exit
    # 130..143) never engages and one large library aborts the whole gated
    # pipeline. Detect that signature and re-exit 137 so the retry kicks in at
    # the next memory tier. Same swallowed-OOM trap we guard in GetOrganelle.
    rc=0
    OMP_NUM_THREADS=${task.cpus} \\
    respect \\
        --debug \\
        --threads ${task.cpus} \\
        --decomp gzip \\
        -d respect_in \\
        -o respect_out > respect_run.log 2>&1 || rc=\$?
    cat respect_run.log

    if [ "\$rc" -ne 0 ]; then
        if grep -qiE "it.s skipped|Number of processes must be at least 1|Killed|MemoryError|bad_alloc|Cannot allocate memory" respect_run.log; then
            echo "[genome_size] RESPECT exit \$rc with a skipped-input/OOM signature -> re-exiting 137 so the GWDG memory-scaling retry engages" >&2
            exit 137
        fi
        echo "[genome_size] RESPECT exit \$rc with no OOM signature -> propagating failure" >&2
        exit \$rc
    fi

    # RESPECT writes estimated-parameters.txt: tab-separated, header row
    # then one row per input. Columns include:
    #   sample, sequencing_depth, genome_length, ...
    # Pull the headline numbers for downstream consumers.
    python3 - <<'PY'
import csv, os, pathlib
out_dir = pathlib.Path('respect_out')
est = out_dir / 'estimated-parameters.txt'
size_v, cov_v = 'NA', 'NA'
with est.open() as fh:
    rows = list(csv.DictReader(fh, delimiter='\\t'))
if rows:
    r = rows[0]
    # RESPECT column names have varied across versions. Try the common ones.
    size_v = r.get('genome_length') or r.get('genome_size') or 'NA'
    cov_v  = r.get('sequencing_depth') or r.get('coverage') or 'NA'

sample = '${meta.id}'
expected = '${expected}'
with open(f'{sample}.genome_size.txt', 'w') as fh:
    fh.write(f'{sample}\\t{size_v}\\t{expected}\\n')
with open(f'{sample}.coverage.txt', 'w') as fh:
    fh.write(f'{sample}\\t{cov_v}\\n')

# MultiQC custom-content TSV. Header magic comments tell MultiQC to render
# this as a table.
with open(f'{sample}_respect_mqc.tsv', 'w') as fh:
    fh.write('# id: respect_genome_size\\n')
    fh.write('# section_name: "RESPECT (genome size + coverage)"\\n')
    fh.write('# description: "Skim-data estimates from RESPECT. expected_size_bp is the samplesheet override (NA if not provided)."\\n')
    fh.write('# format: "tsv"\\n')
    fh.write('# plot_type: "table"\\n')
    fh.write('Sample\\tgenome_size_bp\\tcoverage_x\\texpected_size_bp\\n')
    fh.write(f'{sample}\\t{size_v}\\t{cov_v}\\t{expected}\\n')
PY
    """
}

workflow GENOME_SIZE {
    take:
    reads_ch       // channel: tuple(meta, reads)

    main:
    RESPECT(reads_ch)

    emit:
    size     = RESPECT.out.size
    coverage = RESPECT.out.coverage
    summary  = RESPECT.out.summary
}
