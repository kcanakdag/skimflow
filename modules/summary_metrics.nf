// Per-sample metrics table. Two steps:
//   ASSEMBLY_STATS  - cheap contig stats (n_contigs, total_bp, N50, largest)
//                     per assembly, so the join below doesn't have to stage the
//                     (large) contig FASTAs.
//   SUMMARY_METRICS - join every stage's small summary into one wide CSV, plus
//                     a MultiQC custom-content table. Parsing lives in
//                     bin/collect_metrics.py (stdlib; auto-added to PATH).

process ASSEMBLY_STATS {
    tag "${meta.id}"
    container params.gene_container
    publishDir "${params.outdir}/assembly", mode: 'copy', saveAs: { fn -> "${meta.id}/${fn}" }

    cpus 1
    memory '2 GB'
    time '30m'
    errorStrategy 'ignore'

    input:
    tuple val(meta), path(contigs)

    output:
    tuple val(meta), path("${meta.id}.assembly_stats.tsv"), emit: stats

    script:
    """
    python3 - <<'PY'
lengths = []
seq = 0
with open("${contigs}") as fh:
    for line in fh:
        if line.startswith('>'):
            if seq:
                lengths.append(seq)
            seq = 0
        else:
            seq += len(line.strip())
    if seq:
        lengths.append(seq)

lengths.sort(reverse=True)
total = sum(lengths)
n = len(lengths)
largest = lengths[0] if lengths else 0
half = total / 2.0
acc = 0
n50 = 0
for L in lengths:
    acc += L
    if acc >= half:
        n50 = L
        break

with open("${meta.id}.assembly_stats.tsv", "w") as out:
    out.write("sample\\tn_contigs\\ttotal_bp\\tN50\\tlargest\\n")
    out.write("${meta.id}\\t%d\\t%d\\t%d\\t%d\\n" % (n, total, n50, largest))
PY
    """
}

process SUMMARY_METRICS {
    tag 'summary'
    container params.gene_container
    publishDir "${params.outdir}/summary", mode: 'copy'

    cpus 1
    memory '2 GB'
    time '30m'

    input:
    path qc,      stageAs: 'qc/*'
    path respect, stageAs: 'respect/*'
    path asm,     stageAs: 'assembly/*'
    path mito,    stageAs: 'mito/*'
    path annot,   stageAs: 'annot/*'
    path occ,     stageAs: 'occ/*'
    path busco,   stageAs: 'busco/*'

    output:
    path 'skimflow_metrics.csv',     emit: csv
    path 'skimflow_metrics_mqc.tsv', emit: summary

    script:
    """
    # stageAs only creates a dir when at least one file is present; make the
    # rest so collect_metrics.py always sees all seven inputs (same trick as
    # the REPORT step).
    mkdir -p qc respect assembly mito annot occ busco
    collect_metrics.py \\
        --qc-dir qc \\
        --respect-dir respect \\
        --assembly-dir assembly \\
        --mito-dir mito \\
        --annot-dir annot \\
        --occ-dir occ \\
        --busco-dir busco \\
        --outdir .
    """
}
