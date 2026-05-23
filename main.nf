#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { READ_QC }            from './modules/read_qc.nf'
include { GENOME_SIZE }        from './modules/genome_size.nf'
include { ASSEMBLY }           from './modules/assembly.nf'
include { MITO }               from './modules/mitogenome.nf'
include { MARKER_EXTRACTION }  from './modules/markers.nf'
include { REPORT }             from './modules/report.nf'

// Parse samplesheet.csv into a channel of:
//   tuple(sample_id, [fastq_1, fastq_2], species_id)
// Single-end is detected when fastq_2 is empty.
def parseSamplesheet(csvPath) {
    Channel
        .fromPath(csvPath, checkIfExists: true)
        .splitCsv(header: true, strip: true)
        .map { row ->
            if (!row.sample_id?.trim()) error "samplesheet: missing sample_id in row: ${row}"
            if (!row.fastq_1?.trim())  error "samplesheet: missing fastq_1 for ${row.sample_id}"

            // Relative paths in the CSV are resolved against projectDir.
            def f1     = row.fastq_1.startsWith('/') ? row.fastq_1 : "${projectDir}/${row.fastq_1}"
            def r1     = file(f1, checkIfExists: true)
            def r2_str = row.fastq_2?.trim()
            def f2     = r2_str ? (r2_str.startsWith('/') ? r2_str : "${projectDir}/${r2_str}") : null
            def r2     = f2 ? file(f2, checkIfExists: true) : null
            def reads  = r2 ? [r1, r2] : [r1]
            def size_str = row.expected_genome_size_bp?.trim()
            def expected_size = (size_str && size_str.isLong()) ? size_str.toLong() : null
            def meta = [
                id:            row.sample_id,
                species:       row.species_id ?: 'unknown',
                single_end:    r2 == null,
                expected_size: expected_size,
            ]

            tuple(meta, reads)
        }
}

// Build a one-row channel from --r1/--r2 flags (single-sample shortcut).
// sample_id falls back to R1's basename with a trailing _R1/_1 stripped, so
// Pmisa_R1.fastq.gz -> Pmisa, sample.1.fq.gz -> sample.
def parseDirectReads() {
    def r1 = file(params.r1, checkIfExists: true)
    def r2 = params.r2 ? file(params.r2, checkIfExists: true) : null
    def reads = r2 ? [r1, r2] : [r1]

    def derived_id = r1.simpleName.replaceAll(/[._-]?[Rr]?1$/, '')
    def id         = params.sample_id ?: (derived_id ?: r1.simpleName)
    def size_v     = params.expected_size ? params.expected_size.toString().toLong() : null

    def meta = [
        id:            id,
        species:       params.species ?: id,
        single_end:    r2 == null,
        expected_size: size_v,
    ]
    Channel.of(tuple(meta, reads))
}

workflow {
    if (params.r1) {
        reads_ch = parseDirectReads()
    } else if (params.input) {
        reads_ch = parseSamplesheet(params.input)
    } else {
        error """No input given. Pick one:
  --r1 <R1.fq.gz> [--r2 <R2.fq.gz>]   single-sample shortcut
  --input <samplesheet.csv>           multi-sample (see assets/samplesheet_test.csv)
"""
    }

    READ_QC(reads_ch)
    ASSEMBLY(READ_QC.out.reads)
    MITO(READ_QC.out.reads)
    MARKER_EXTRACTION(ASSEMBLY.out.contigs)

    if (params.skip_respect) {
        genome_size_ch = Channel.empty().collect().ifEmpty([])
    } else {
        GENOME_SIZE(READ_QC.out.reads)
        genome_size_ch = GENOME_SIZE.out.summary.map { meta, f -> f }.collect()
    }

    REPORT(
        READ_QC.out.json.map              { meta, f -> f }.collect(),
        genome_size_ch,
        MARKER_EXTRACTION.out.summary.map { meta, f -> f }.collect().ifEmpty([]),
    )
}
