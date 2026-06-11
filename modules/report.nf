// Final per-run report. MultiQC scrapes its supported upstream outputs
// (fastp.json, RESPECT summary, BUSCO short_summary, kraken2 reports, and
// skimflow mitogenome custom-content summaries) and rolls them into one HTML.

process REPORT {
    tag 'multiqc'
    container 'quay.io/biocontainers/multiqc:1.34--pyhdfd78af_0'
    publishDir "${params.outdir}/report", mode: 'copy'

    cpus 1
    memory '2 GB'

    input:
    path qc_files,          stageAs: 'qc/*'
    path genome_size_files, stageAs: 'genome_size/*'
    path busco_summaries,   stageAs: 'busco/*'
    path kraken_reports,    stageAs: 'kraken/*'
    path mitogenome_summaries, stageAs: 'mitogenome/*'

    output:
    path 'multiqc_report.html',          emit: html
    path 'multiqc_report_data',          emit: data

    script:
    """
    # Nextflow's stageAs only creates a directory when there's at least one
    # input file; if a step was skipped (e.g. --skip_respect) or all its tasks
    # failed-ignored (mito/markers on tiny test data), the dir is absent and
    # multiqc bails out. Create empty placeholders to keep the CLI uniform.
    mkdir -p qc genome_size busco kraken mitogenome
cat > multiqc_config.yaml <<'EOF'
custom_content:
  order:
    - getorganelle_mitogenome
    - mitoz_annotation
    - mitos2_annotation
EOF

    multiqc \\
        --force \\
        --config multiqc_config.yaml \\
        --filename multiqc_report.html \\
        --title 'skimflow' \\
        qc genome_size busco kraken mitogenome
    """
}
