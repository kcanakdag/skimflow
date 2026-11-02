// De novo assembly of short reads with MEGAHIT.
// MEGAHIT is fast and memory-frugal; --presets meta-sensitive matches the
// host + microbial mixture typical of genome-skim data. The output dir must
// NOT pre-exist, which a fresh Nextflow task dir guarantees.

process ASSEMBLY {
    tag "${meta.id}"
    container 'quay.io/biocontainers/megahit:1.2.9--haf24da9_8'
    publishDir "${params.outdir}/assembly", mode: 'copy', saveAs: { fn -> "${meta.id}/${fn}" }

    cpus 16
    // Scale with task.attempt so the gwdg retry on exit 130..143 (OOM or
    // walltime kill) gets more room instead of repeating the same request.
    memory { 32.GB * task.attempt }
    // MEGAHIT's meta-sensitive preset runs the full k-ladder (k=21..141);
    // larger genomes (~1 Gbp skim) exceed 4 h and get wall-clock-killed
    // (exit 140). scc-cpu allows up to 48 h; 24 h gives ample headroom and
    // a retry is capped at 47 h to stay under that limit.
    time { [24.h * task.attempt, 47.h].min() }

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("${meta.id}.contigs.fasta"), emit: contigs
    tuple val(meta), path("${meta.id}.megahit.log"),   emit: logfile

    script:
    def input_args = meta.single_end ? "-r ${reads[0]}" : "-1 ${reads[0]} -2 ${reads[1]}"
    """
    megahit \\
        ${input_args} \\
        --presets ${params.megahit_preset} \\
        -t ${task.cpus} \\
        -m ${(task.memory.toBytes() * 0.9) as long} \\
        -o megahit_out

    cp megahit_out/final.contigs.fa ${meta.id}.contigs.fasta
    cp megahit_out/log              ${meta.id}.megahit.log
    """
}
