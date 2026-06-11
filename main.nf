#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { READ_QC }               from './modules/read_qc.nf'
include { DECONTAM }              from './modules/decontam.nf'
include { GENOME_SIZE }           from './modules/genome_size.nf'
include { ASSEMBLY }              from './modules/assembly.nf'
include { LR_QC }                 from './modules/long_read_qc.nf'
include { FLYE }                  from './modules/flye.nf'
include { MITO }                  from './modules/mitogenome.nf'
include { MITOGENOME_ANNOTATION } from './modules/mitogenome_annotation.nf'
include { MARKER_EXTRACTION }     from './modules/markers.nf'
include { REPORT }                from './modules/report.nf'

// Allowed long-read technologies (used to validate meta.lr_type at parse time).
// Declared as a function (not a top-level statement) so the strict Nextflow
// syntax parser accepts it alongside the other declarations.
def LR_TYPES() { ['nanopore', 'pacbio'] }

// Parse samplesheet.csv into a channel of tuple(meta, reads, lr).
// Exactly one of reads (non-empty list => short) / lr (path => long) is
// populated; meta.read_type records which. Two new optional columns are
// backward compatible via splitCsv(header:true): long_reads, lr_type.
def parseSamplesheet(csvPath) {
    Channel
        .fromPath(csvPath, checkIfExists: true)
        .splitCsv(header: true, strip: true)
        .map { row ->
            if (!row.sample_id?.trim()) error "samplesheet: missing sample_id in row: ${row}"

            def has_short = row.fastq_1?.trim() as boolean
            def lr_str    = row.long_reads?.trim()
            def has_long  = lr_str as boolean

            if (has_short && has_long)
                error "samplesheet: ${row.sample_id} has both fastq_1 and long_reads; a sample is short OR long, not both"
            if (!has_short && !has_long)
                error "samplesheet: ${row.sample_id} has neither fastq_1 nor long_reads"

            def size_str      = row.expected_genome_size_bp?.trim()
            def expected_size = (size_str && size_str.isLong()) ? size_str.toLong() : null

            if (has_long) {
                def lf  = lr_str.startsWith('/') ? lr_str : "${projectDir}/${lr_str}"
                def lr  = file(lf, checkIfExists: true)
                def lrt = (row.lr_type?.trim() ?: params.lr_type)
                if (!(lrt in LR_TYPES()))
                    error "samplesheet: ${row.sample_id} lr_type '${lrt}' must be one of ${LR_TYPES()}"
                def meta = [
                    id:            row.sample_id,
                    species:       row.species_id ?: 'unknown',
                    single_end:    true,
                    expected_size: expected_size,
                    read_type:     'long',
                    lr_type:       lrt,
                ]
                return tuple(meta, [], lr)
            }

            // Short-read sample (the original short-read path).
            def f1     = row.fastq_1.startsWith('/') ? row.fastq_1 : "${projectDir}/${row.fastq_1}"
            def r1     = file(f1, checkIfExists: true)
            def r2_str = row.fastq_2?.trim()
            def f2     = r2_str ? (r2_str.startsWith('/') ? r2_str : "${projectDir}/${r2_str}") : null
            def r2     = f2 ? file(f2, checkIfExists: true) : null
            def reads  = r2 ? [r1, r2] : [r1]
            def meta = [
                id:            row.sample_id,
                species:       row.species_id ?: 'unknown',
                single_end:    r2 == null,
                expected_size: expected_size,
                read_type:     'short',
                lr_type:       null,
            ]
            tuple(meta, reads, null)
        }
}

// Build a one-row channel from --r1/--r2 flags (short single-sample shortcut).
// sample_id falls back to R1's basename with a trailing _R1/_1 stripped.
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
        read_type:     'short',
        lr_type:       null,
    ]
    Channel.of(tuple(meta, reads, null))
}

// Build a one-row channel from --long_reads (long single-sample shortcut).
def parseDirectLong() {
    def lr  = file(params.long_reads, checkIfExists: true)
    def lrt = params.lr_type
    if (!(lrt in LR_TYPES()))
        error "--lr_type '${lrt}' must be one of ${LR_TYPES()}"

    def id     = params.sample_id ?: lr.simpleName
    def size_v = params.expected_size ? params.expected_size.toString().toLong() : null

    def meta = [
        id:            id,
        species:       params.species ?: id,
        single_end:    true,
        expected_size: size_v,
        read_type:     'long',
        lr_type:       lrt,
    ]
    Channel.of(tuple(meta, [], lr))
}

workflow {
    // Mutual exclusion must be an explicit guard: an if/else if ladder would
    // silently pick the first branch and never reach a "both given" error.
    if (params.long_reads && (params.r1 || params.input))
        error "--long_reads is long-only direct mode; do not combine with --r1 or --input"

    if (params.r1) {
        parsed_ch = parseDirectReads()       // short
    } else if (params.long_reads) {
        parsed_ch = parseDirectLong()        // long
    } else if (params.input) {
        parsed_ch = parseSamplesheet(params.input)   // either / both (per row)
    } else {
        error """No input given. Pick one:
  --r1 <R1.fq.gz> [--r2 <R2.fq.gz>]   single-sample short-read shortcut
  --long_reads <reads.fq.gz>          single-sample long-read shortcut
  --input <samplesheet.csv>           multi-sample (see assets/samplesheet_test.csv)
"""
    }
    // Split into the two disjoint tracks. No join/remainder: samples never
    // belong to both, so the two contigs streams simply mix before BUSCO.
    short_in = parsed_ch.filter { meta, reads, lr -> meta.read_type == 'short' }
                        .map    { meta, reads, lr -> tuple(meta, reads) }
    long_in  = parsed_ch.filter { meta, reads, lr -> meta.read_type == 'long' }
                        .map    { meta, reads, lr -> tuple(meta, lr) }

    // --- Short-read track ---
    READ_QC(short_in)

    if (params.kraken2_db) {
        kraken2_db_ch = Channel.value(file(params.kraken2_db, checkIfExists: true))
        DECONTAM(READ_QC.out.reads, kraken2_db_ch)
        clean_reads      = DECONTAM.out.reads
        kraken_report_ch = DECONTAM.out.report.map { meta, f -> f }.collect().ifEmpty([])
    } else {
        clean_reads      = READ_QC.out.reads
        kraken_report_ch = Channel.empty().collect().ifEmpty([])
    }

    ASSEMBLY(clean_reads)        // MEGAHIT
    MITO(clean_reads)
    MITOGENOME_ANNOTATION(MITO.out.fasta)
    mitogenome_summary_ch = MITO.out.summary
        .mix(MITOGENOME_ANNOTATION.out.summary)
        .map { meta, f -> f }
        .collect()
        .ifEmpty([])

    if (params.skip_respect) {
        genome_size_ch = Channel.empty().collect().ifEmpty([])
    } else {
        GENOME_SIZE(clean_reads)
        genome_size_ch = GENOME_SIZE.out.summary.map { meta, f -> f }.collect().ifEmpty([])
    }

    // --- Long-read track ---
    LR_QC(long_in)
    FLYE(LR_QC.out.reads)

    // --- Markers over both assemblies (mix is where the tracks rejoin) ---
    MARKER_EXTRACTION(ASSEMBLY.out.contigs.mix(FLYE.out.contigs))

    REPORT(
        READ_QC.out.json.map              { meta, f -> f }.collect().ifEmpty([]),
        genome_size_ch,
        MARKER_EXTRACTION.out.summary.map { meta, f -> f }.collect().ifEmpty([]),
        kraken_report_ch,
        mitogenome_summary_ch,
    )
}
