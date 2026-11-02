#!/usr/bin/env nextflow

include { READ_QC }               from './modules/read_qc.nf'
include { DECONTAM }              from './modules/decontam.nf'
include { GENOME_SIZE }           from './modules/genome_size.nf'
include { ASSEMBLY }              from './modules/assembly.nf'
include { LR_QC }                 from './modules/long_read_qc.nf'
include { FLYE }                  from './modules/flye.nf'
include { MITO }                  from './modules/mitogenome.nf'
include { MITOGENOME_ANNOTATION } from './modules/mitogenome_annotation.nf'
include { MITO_GENES }            from './modules/mito_genes.nf'
include { MARKER_EXTRACTION }     from './modules/markers.nf'
include { PHYLOGENY }             from './modules/phylogeny.nf'
include { ASSEMBLY_STATS }        from './modules/summary_metrics.nf'
include { SUMMARY_METRICS }       from './modules/summary_metrics.nf'
include { REPORT }                from './modules/report.nf'

// Allowed long-read technologies (used to validate meta.lr_type at parse time).
// Declared as a function (not a top-level statement) so the strict Nextflow
// syntax parser accepts it alongside the other declarations.
def LR_TYPES() { ['nanopore', 'pacbio'] }

// Sample ids name every per-sample output file and results/ subdirectory and
// are passed to tool command lines, so restrict them to a file-name-safe set.
def SAMPLE_ID_RE() { /^[A-Za-z0-9_-]+$/ }

def checkSampleId(id, where) {
    if (!(id ==~ SAMPLE_ID_RE()))
        error "${where}: sample_id '${id}' may only contain letters, digits, '_' and '-' (it names output files; replace dots or spaces with underscores)"
    id
}

// Resolve a samplesheet path: absolute paths and URLs (s3://, https://, ...)
// are used as-is, a leading ~/ expands to $HOME, and any other relative path
// is resolved against projectDir (the documented rule).
def resolveSheetPath(String p) {
    if (p.startsWith('/') || p.contains('://')) return p
    if (p.startsWith('~/')) return System.getenv('HOME') + p.substring(1)
    "${projectDir}/${p}"
}

// Parse an expected genome size in bp. Blank means "not given" (null). Whole
// numbers in plain or scientific notation (1000000000, 1e9) are accepted;
// anything else (1.5, 1,000, text) is a clear error instead of a silent NA.
def parseExpectedSize(v, where) {
    if (v == null || v.toString().trim() == '') return null
    try {
        return new BigDecimal(v.toString().trim()).toBigIntegerExact().longValueExact()
    } catch (Exception e) {
        error "${where}: expected_genome_size_bp must be a whole number of bp, got '${v}'"
    }
}

// Parse samplesheet.csv into a channel of tuple(meta, reads, lr).
// Exactly one of reads (non-empty list => short) / lr (path => long) is
// populated; meta.read_type records which. Two new optional columns are
// backward compatible via splitCsv(header:true): long_reads, lr_type.
def parseSamplesheet(csvPath) {
    // Validate ids eagerly so a bad sheet fails before any task runs, not at
    // the cross-sample steps after all per-sample compute.
    def ids  = file(csvPath, checkIfExists: true).splitCsv(header: true, strip: true)*.sample_id*.trim().findAll { it }
    def dups = ids.countBy { it }.findAll { k, v -> v > 1 }.keySet()
    if (dups) error "samplesheet: duplicate sample_id(s): ${dups.join(', ')}"
    ids.each { checkSampleId(it, 'samplesheet') }

    Channel
        .fromPath(csvPath, checkIfExists: true)
        .splitCsv(header: true, strip: true)
        .map { row ->
            if (!row.sample_id?.trim()) error "samplesheet: missing sample_id in row: ${row}"
            def sid = row.sample_id.trim()

            def has_short = row.fastq_1?.trim() as boolean
            def lr_str    = row.long_reads?.trim()
            def has_long  = lr_str as boolean

            if (has_short && has_long)
                error "samplesheet: ${row.sample_id} has both fastq_1 and long_reads; a sample is short OR long, not both"
            if (!has_short && !has_long)
                error "samplesheet: ${row.sample_id} has neither fastq_1 nor long_reads"

            def size_str      = row.expected_genome_size_bp?.trim()
            def expected_size = parseExpectedSize(size_str, "samplesheet: ${sid}")

            if (has_long) {
                def lf  = resolveSheetPath(lr_str)
                def lr  = file(lf, checkIfExists: true)
                def lrt = (row.lr_type?.trim() ?: params.lr_type)
                if (!(lrt in LR_TYPES()))
                    error "samplesheet: ${row.sample_id} lr_type '${lrt}' must be one of ${LR_TYPES()}"
                def meta = [
                    id:            sid,
                    species:       row.species_id?.trim() ?: sid,
                    single_end:    true,
                    expected_size: expected_size,
                    read_type:     'long',
                    lr_type:       lrt,
                ]
                return tuple(meta, [], lr)
            }

            // Short-read sample (the original short-read path).
            def f1     = resolveSheetPath(row.fastq_1.trim())
            def r1     = file(f1, checkIfExists: true)
            def r2_str = row.fastq_2?.trim()
            def f2     = r2_str ? resolveSheetPath(r2_str) : null
            def r2     = f2 ? file(f2, checkIfExists: true) : null
            def reads  = r2 ? [r1, r2] : [r1]
            def meta = [
                id:            sid,
                species:       row.species_id?.trim() ?: sid,
                single_end:    r2 == null,
                expected_size: expected_size,
                read_type:     'short',
                lr_type:       null,
            ]
            tuple(meta, reads, null)
        }
}

// Build a one-row channel from --r1/--r2 flags (short single-sample shortcut).
// sample_id falls back to R1's basename with a trailing _R1/_1 (optionally
// followed by Illumina's _001) stripped; a bare trailing digit such as in
// sample01 or S11 is kept.
def parseDirectReads() {
    def r1 = file(params.r1, checkIfExists: true)
    def r2 = params.r2 ? file(params.r2, checkIfExists: true) : null
    def reads = r2 ? [r1, r2] : [r1]

    def derived_id = r1.simpleName.replaceAll(/([._-]?[Rr]1|[._-]1)(_001)?$/, '')
    def id         = checkSampleId((params.sample_id ?: (derived_id ?: r1.simpleName)).toString().trim(), '--r1')
    def size_v     = parseExpectedSize(params.expected_size ?: null, '--expected_size')

    def meta = [
        id:            id,
        species:       params.species?.toString()?.trim() ?: id,
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

    def id     = checkSampleId((params.sample_id ?: lr.simpleName).toString().trim(), '--long_reads')
    def size_v = parseExpectedSize(params.expected_size ?: null, '--expected_size')

    def meta = [
        id:            id,
        species:       params.species?.toString()?.trim() ?: id,
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
    if (params.r1 && params.input)
        error "--r1/--r2 is single-sample direct mode; do not combine with --input"
    if (params.r2 && !params.r1)
        error "--r2 requires --r1"
    if (!(params.mitogenome_topology in ['auto', 'linear', 'circular']))
        error "mitogenome_topology must be one of: auto, linear, circular (got '${params.mitogenome_topology}')"

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
        kraken_report_ch = Channel.value([])
    }

    ASSEMBLY(clean_reads)        // MEGAHIT
    MITO(clean_reads)
    MITOGENOME_ANNOTATION(MITO.out.fasta)
    mitogenome_summary_ch = MITO.out.summary
        .mix(MITOGENOME_ANNOTATION.out.summary)
        .map { meta, f -> f }
        .collect()
        .ifEmpty([])

    MITO_GENES(MITOGENOME_ANNOTATION.out.mitos2_dir)
    gene_summary_ch = MITO_GENES.out.summary.collect().ifEmpty([])

    if (params.skip_respect) {
        genome_size_ch = Channel.value([])
    } else {
        GENOME_SIZE(clean_reads)
        genome_size_ch = GENOME_SIZE.out.summary.map { meta, f -> f }.collect().ifEmpty([])
    }

    // --- Long-read track ---
    LR_QC(long_in)
    FLYE(LR_QC.out.reads)

    // --- Markers over both assemblies (mix is where the tracks rejoin) ---
    contigs_ch = ASSEMBLY.out.contigs.mix(FLYE.out.contigs)
    MARKER_EXTRACTION(contigs_ch)
    ASSEMBLY_STATS(contigs_ch)

    // --- Phylogeny over the harvested per-gene FASTAs (nt + aa) ---
    if (params.skip_phylo) {
        phylo_summary_ch = Channel.empty()
        phylo_tree_ch    = Channel.empty()
    } else {
        PHYLOGENY(MITO_GENES.out.nt, MITO_GENES.out.aa)
        phylo_summary_ch = PHYLOGENY.out.summary
        phylo_tree_ch    = PHYLOGENY.out.treeplot
    }

    // --- Per-sample metrics table joining every stage ---
    // Samplesheet metadata (id, species, expected size) for both tracks, so
    // the table still has species / expected size when a stage is skipped.
    meta_tsv_ch = parsed_ch
        .map { meta, reads, lr -> "${meta.id}\t${meta.species ?: ''}\t${meta.expected_size ?: ''}" }
        .collectFile(name: 'samples_meta.tsv', newLine: true)
        .ifEmpty([])
    // QC json and BUSCO summaries feed both SUMMARY_METRICS and REPORT.
    qc_json_ch = READ_QC.out.json.map { meta, f -> f }.collect().ifEmpty([])
    busco_ch   = MARKER_EXTRACTION.out.summary.map { meta, f -> f }.collect().ifEmpty([])

    SUMMARY_METRICS(
        qc_json_ch,
        genome_size_ch,
        ASSEMBLY_STATS.out.stats.map            { meta, f -> f }.collect().ifEmpty([]),
        MITO.out.summary.map                    { meta, f -> f }.collect().ifEmpty([]),
        MITOGENOME_ANNOTATION.out.summary.map   { meta, f -> f }.collect().ifEmpty([]),
        MITO_GENES.out.matrix.collect().ifEmpty([]),
        busco_ch,
        meta_tsv_ch,
    )

    // Phylogeny tables, the rendered tree SVGs, and the metrics table all feed
    // MultiQC's extra dir as custom content.
    extra_summaries_ch = phylo_summary_ch
        .mix(phylo_tree_ch)
        .mix(SUMMARY_METRICS.out.summary)
        .collect()
        .ifEmpty([])

    REPORT(
        qc_json_ch,
        genome_size_ch,
        busco_ch,
        kraken_report_ch,
        mitogenome_summary_ch,
        gene_summary_ch,
        extra_summaries_ch,
    )
}
