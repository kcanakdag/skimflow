// Read QC: fastp for adapter trimming + quality filtering.

process READ_QC {
    tag "${meta.id}"
    container 'quay.io/biocontainers/fastp:1.3.3--h43da1c4_0'
    // fastp's trimmed FASTQs are large intermediates (tens of GB on real
    // libraries) that downstream steps consume from work/, not results/. Publish
    // only the QC reports (json/html) by default; a saveAs returning null skips
    // the reads. Opt back in with --save_trimmed_reads (e.g. to keep the trimmed
    // set as a deliverable). Publishing them by default overran the GWDG home
    // quota and crashed the run mid-flight.
    publishDir "${params.outdir}/read_qc", mode: 'copy',
        saveAs: { fn -> (fn.endsWith('.fastq.gz') && !params.save_trimmed_reads) ? null : fn }

    cpus 2
    memory '2 GB'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("${meta.id}.trimmed_R*.fastq.gz"), emit: reads
    tuple val(meta), path("${meta.id}.fastp.json"),          emit: json
    tuple val(meta), path("${meta.id}.fastp.html"),          emit: html

    script:
    if (meta.single_end) {
        """
        fastp \\
            --in1 ${reads[0]} \\
            --out1 ${meta.id}.trimmed_R1.fastq.gz \\
            --json ${meta.id}.fastp.json \\
            --html ${meta.id}.fastp.html \\
            --thread ${task.cpus}
        """
    } else {
        """
        fastp \\
            --in1 ${reads[0]} --in2 ${reads[1]} \\
            --out1 ${meta.id}.trimmed_R1.fastq.gz \\
            --out2 ${meta.id}.trimmed_R2.fastq.gz \\
            --json ${meta.id}.fastp.json \\
            --html ${meta.id}.fastp.html \\
            --thread ${task.cpus}
        """
    }
}
