#!/usr/bin/env python3
"""Join every stage's per-sample summary into one wide metrics table.

Reads the small summary files each pipeline stage already publishes (staged
here into one subdir per stage) and emits, per sample:
  results/summary/skimflow_metrics.csv        one row per sample, all columns
  results/summary/skimflow_metrics_mqc.tsv    curated subset for MultiQC

Every parser is defensive: a stage that was skipped or failed-soft simply
leaves its columns blank rather than aborting the join. Standard library only.
"""
import argparse
import csv
import glob
import json
import os
import re
import sys

# Column order of the wide CSV.
COLUMNS = [
    'sample', 'species',
    'reads_in', 'reads_out', 'pct_q30',
    'genome_size_bp', 'coverage_x', 'expected_size_bp',
    'assembly_contigs', 'assembly_total_bp', 'assembly_N50',
    'mito_length_bp', 'mito_status',
    'mitoz_features', 'mitos2_features',
    'genes_recovered',
    'busco_C', 'busco_S', 'busco_D', 'busco_F', 'busco_M',
]
GENE_STEMS = ['COX1', 'COX2', 'COX3', 'COB', 'NAD1', 'NAD2', 'NAD3', 'NAD4',
              'NAD4L', 'NAD5', 'NAD6', 'ATP6', 'ATP8', '16S', '12S']


def _rows(d):
    """Yield rows (as dicts) from a headered TSV, skipping MultiQC # comments."""
    lines = [ln for ln in d.splitlines() if ln.strip() and not ln.startswith('#')]
    if not lines:
        return
    reader = csv.DictReader(lines, delimiter='\t')
    for row in reader:
        yield row


def _read(path):
    try:
        with open(path) as fh:
            return fh.read()
    except OSError:
        return ''


def rec(table, sample):
    return table.setdefault(sample, {'sample': sample})


def parse_qc(qc_dir, table):
    for path in sorted(glob.glob(os.path.join(qc_dir, '*.fastp.json'))):
        sample = os.path.basename(path)[:-len('.fastp.json')]
        try:
            data = json.loads(_read(path))
            r = rec(table, sample)
            before = data.get('summary', {}).get('before_filtering', {})
            after = data.get('summary', {}).get('after_filtering', {})
            r['reads_in'] = before.get('total_reads', '')
            r['reads_out'] = after.get('total_reads', '')
            q30 = after.get('q30_rate', '')
            r['pct_q30'] = round(q30 * 100, 2) if isinstance(q30, (int, float)) else ''
        except (ValueError, KeyError, TypeError):
            continue


def parse_respect(respect_dir, table):
    for path in sorted(glob.glob(os.path.join(respect_dir, '*_respect_mqc.tsv'))):
        for row in _rows(_read(path)):
            sample = row.get('Sample') or os.path.basename(path)[:-len('_respect_mqc.tsv')]
            r = rec(table, sample)
            r['genome_size_bp'] = row.get('genome_size_bp', '')
            r['coverage_x'] = row.get('coverage_x', '')
            r['expected_size_bp'] = row.get('expected_size_bp', '')


def parse_assembly(asm_dir, table):
    for path in sorted(glob.glob(os.path.join(asm_dir, '*.assembly_stats.tsv'))):
        for row in _rows(_read(path)):
            sample = row.get('sample') or os.path.basename(path)[:-len('.assembly_stats.tsv')]
            r = rec(table, sample)
            r['assembly_contigs'] = row.get('n_contigs', '')
            r['assembly_total_bp'] = row.get('total_bp', '')
            r['assembly_N50'] = row.get('N50', '')


def parse_mito(mito_dir, table):
    for path in sorted(glob.glob(os.path.join(mito_dir, '*.getorganelle_mqc.tsv'))):
        for row in _rows(_read(path)):
            sample = row.get('Sample')
            if not sample:
                continue
            r = rec(table, sample)
            r['mito_length_bp'] = row.get('Longest_bp', '')
            r['mito_status'] = row.get('Status', '')


def parse_annot(annot_dir, table):
    for path in sorted(glob.glob(os.path.join(annot_dir, '*.mitoz_mqc.tsv'))):
        for row in _rows(_read(path)):
            if row.get('Sample'):
                rec(table, row['Sample'])['mitoz_features'] = row.get('Features', '')
    for path in sorted(glob.glob(os.path.join(annot_dir, '*.mitos2_mqc.tsv'))):
        for row in _rows(_read(path)):
            if row.get('Sample'):
                rec(table, row['Sample'])['mitos2_features'] = row.get('Features', '')


def parse_occupancy(occ_dir, table):
    for path in sorted(glob.glob(os.path.join(occ_dir, 'occupancy.tsv'))):
        for row in _rows(_read(path)):
            sample = row.get('sample')
            if not sample:
                continue
            r = rec(table, sample)
            if row.get('species'):
                r['species'] = row['species']
            present = 0
            for g in GENE_STEMS:
                v = (row.get(g) or '').strip()
                if v and v != '0':
                    present += 1
            r['genes_recovered'] = present


_BUSCO_RE = re.compile(
    r'C:([\d.]+)%\[S:([\d.]+)%,D:([\d.]+)%\],F:([\d.]+)%,M:([\d.]+)%')
_BUSCO_FILE_RE = re.compile(r'for file\s+(\S+)')


def parse_busco(busco_dir, table):
    for path in sorted(glob.glob(os.path.join(busco_dir, 'short_summary*'))):
        text = _read(path)
        m = _BUSCO_RE.search(text)
        if not m:
            continue
        sample = None
        fm = _BUSCO_FILE_RE.search(text)
        if fm:
            sample = os.path.basename(fm.group(1)).split('.', 1)[0]
        if not sample:
            # Fall back to the run name embedded in the filename: ...<id>_busco.txt
            base = os.path.basename(path)
            rn = re.search(r'\.([^.]+)_busco\.txt$', base)
            sample = rn.group(1) if rn else base
        r = rec(table, sample)
        r['busco_C'], r['busco_S'], r['busco_D'], r['busco_F'], r['busco_M'] = m.groups()


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument('--qc-dir', default='qc')
    ap.add_argument('--respect-dir', default='respect')
    ap.add_argument('--assembly-dir', default='assembly')
    ap.add_argument('--mito-dir', default='mito')
    ap.add_argument('--annot-dir', default='annot')
    ap.add_argument('--occ-dir', default='occ')
    ap.add_argument('--busco-dir', default='busco')
    ap.add_argument('--outdir', required=True)
    args = ap.parse_args(argv)

    table = {}
    parse_qc(args.qc_dir, table)
    parse_respect(args.respect_dir, table)
    parse_assembly(args.assembly_dir, table)
    parse_mito(args.mito_dir, table)
    parse_annot(args.annot_dir, table)
    parse_occupancy(args.occ_dir, table)
    parse_busco(args.busco_dir, table)

    os.makedirs(args.outdir, exist_ok=True)
    samples = sorted(table)

    csv_path = os.path.join(args.outdir, 'skimflow_metrics.csv')
    with open(csv_path, 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=COLUMNS, extrasaction='ignore')
        w.writeheader()
        for s in samples:
            w.writerow(table[s])

    # Curated MultiQC table (custom content).
    mqc_cols = ['reads_out', 'coverage_x', 'mito_length_bp', 'genes_recovered', 'busco_C']
    mqc_path = os.path.join(args.outdir, 'skimflow_metrics_mqc.tsv')
    with open(mqc_path, 'w') as fh:
        fh.write('# id: "skimflow_metrics"\n')
        fh.write('# section_name: "skimflow per-sample metrics"\n')
        fh.write('# description: "One row per sample joining QC, RESPECT, assembly, mitogenome, gene occupancy and BUSCO. Blank = stage skipped or produced no result."\n')
        fh.write('# plot_type: "table"\n')
        fh.write('Sample\t' + '\t'.join(mqc_cols) + '\n')
        for s in samples:
            r = table[s]
            fh.write(s + '\t' + '\t'.join(str(r.get(c, '')) for c in mqc_cols) + '\n')

    sys.stderr.write('[metrics] wrote %d samples to %s\n' % (len(samples), csv_path))
    return 0


if __name__ == '__main__':
    sys.exit(main())
