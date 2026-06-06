// Optional decontamination of short reads with kraken2 ("keep unclassified").
// Reads matching the kraken2 DB (contaminants) are dropped; the unclassified
// reads (the target organism, absent from a PlusPFP-style DB) are written by
// --unclassified-out and become the clean set feeding MEGAHIT, RESPECT, and
// GetOrganelle. Enabled only when --kraken2_db is set (wired in main.nf).
//
// kraken2 loads the whole DB into RAM by default (~17-20 GB for the 16 GB
// index); on a memory-tight node add --memory-mapping to mmap it off disk
// (slower). kraken2 auto-detects gzipped input; if a build ever rejects it,
// add --gzip-compressed. Output FASTQ is uncompressed, so we gzip it to keep
// downstream channels consistent with fastp's .fastq.gz.

process DECONTAM {
    tag "${meta.id}"
    container 'quay.io/biocontainers/kraken2:2.1.6--pl5321h077b44d_0'
    publishDir "${params.outdir}/decontam", mode: 'copy', saveAs: { fn -> "${meta.id}/${fn}" }

    cpus 16
    memory '32 GB'
    time '4h'

    errorStrategy 'terminate'

    input:
    tuple val(meta), path(reads)
    path db

    output:
    tuple val(meta), path("${meta.id}_decontam*.fastq.gz"), emit: reads
    tuple val(meta), path("${meta.id}.kraken.report"),      emit: report

    script:
    if (meta.single_end) {
        """
        kraken2 --db ${db} --threads ${task.cpus} \\
                --unclassified-out ${meta.id}_decontam.fastq \\
                --output ${meta.id}.kraken.out \\
                --report ${meta.id}.kraken.report \\
                ${reads[0]}
        gzip ${meta.id}_decontam.fastq
        """
    } else {
        """
        # kraken2 replaces the # in --unclassified-out with _1 / _2 for paired reads
        kraken2 --db ${db} --threads ${task.cpus} --paired \\
                --unclassified-out ${meta.id}_decontam_R#.fastq \\
                --output ${meta.id}.kraken.out \\
                --report ${meta.id}.kraken.report \\
                ${reads[0]} ${reads[1]}
        gzip ${meta.id}_decontam_R_1.fastq ${meta.id}_decontam_R_2.fastq
        """
    }
}
