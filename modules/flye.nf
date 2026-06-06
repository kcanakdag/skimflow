// Long-read assembly with Flye. Runs only for long-read samples (consumes the
// Filtlong-filtered reads). The read-type flag derives from meta.lr_type;
// --nano-hq is the recommended default for current ONT chemistry (Guppy5+ /
// R10), with params.flye_mode as an explicit override for legacy --nano-raw or
// PacBio HiFi (--pacbio-hifi). Flye >=2.8 does not require --genome-size, so it
// is omitted for skim-sized targets. The contigs emit shape matches MEGAHIT so
// the two assemblers mix into one BUSCO step.

process FLYE {
    tag "${meta.id}"
    container 'quay.io/biocontainers/flye:2.9.6--py313h7fbb527_1'
    publishDir "${params.outdir}/assembly", mode: 'copy', saveAs: { fn -> "${meta.id}/${fn}" }

    cpus 16
    memory '32 GB'
    time '8h'

    // Mirrors the existing assembly/mito tolerance for thin test data. On real
    // long-read data this should succeed; long-read runs are deliberate, so
    // switching to 'terminate' later is defensible.
    errorStrategy 'ignore'

    input:
    tuple val(meta), path(long_reads)

    output:
    tuple val(meta), path("${meta.id}.contigs.fasta"),      emit: contigs
    tuple val(meta), path("${meta.id}.flye_info.txt"),      emit: info
    tuple val(meta), path("${meta.id}.assembly_graph.gfa"), emit: graph, optional: true

    script:
    def mode = params.flye_mode ?: (meta.lr_type == 'pacbio' ? '--pacbio-raw' : '--nano-hq')
    """
    flye ${mode} ${long_reads} \\
         --out-dir flye_out \\
         --threads ${task.cpus}

    cp flye_out/assembly.fasta    ${meta.id}.contigs.fasta
    cp flye_out/assembly_info.txt ${meta.id}.flye_info.txt
    [ -f flye_out/assembly_graph.gfa ] && cp flye_out/assembly_graph.gfa ${meta.id}.assembly_graph.gfa || true
    """
}
