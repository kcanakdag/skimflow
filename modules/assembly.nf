// De novo assembly with SPAdes.
// `--only-assembler` skips read error correction (we trust fastp upstream)
// and `--cov-cutoff auto` drops low-coverage tip contigs typical of skim data.

process ASSEMBLY {
    tag "${meta.id}"
    container 'quay.io/biocontainers/spades:4.2.0--h8d6e82b_2'
    publishDir "${params.outdir}/assembly", mode: 'copy', saveAs: { fn -> "${meta.id}/${fn}" }

    cpus 16
    memory '32 GB'
    time '4h'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("${meta.id}.contigs.fasta"),         emit: contigs
    tuple val(meta), path("${meta.id}.scaffolds.fasta"),       emit: scaffolds, optional: true
    tuple val(meta), path("${meta.id}.assembly_graph.fastg"),  emit: graph,     optional: true
    tuple val(meta), path("${meta.id}.spades.log"),            emit: log

    script:
    def memGb = (task.memory.toGiga() as int) ?: 8
    def input_args = meta.single_end ? "-s ${reads[0]}" : "-1 ${reads[0]} -2 ${reads[1]}"
    """
    spades.py \\
        ${input_args} \\
        -o spades_out \\
        --only-assembler \\
        --cov-cutoff auto \\
        --threads ${task.cpus} \\
        --memory ${memGb}

    cp spades_out/contigs.fasta        ${meta.id}.contigs.fasta
    cp spades_out/spades.log           ${meta.id}.spades.log
    [ -f spades_out/scaffolds.fasta ]      && cp spades_out/scaffolds.fasta      ${meta.id}.scaffolds.fasta      || true
    [ -f spades_out/assembly_graph.fastg ] && cp spades_out/assembly_graph.fastg ${meta.id}.assembly_graph.fastg || true
    """
}
