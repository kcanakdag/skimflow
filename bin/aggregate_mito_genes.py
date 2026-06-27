#!/usr/bin/env python3
"""Aggregate per-sample harvested mito genes into per-gene multi-FASTAs and a
gene-occupancy matrix. Standard library only."""
import argparse
import os
import sys

GENE_STEMS = ['COX1', 'COX2', 'COX3', 'COB', 'NAD1', 'NAD2', 'NAD3', 'NAD4',
              'NAD4L', 'NAD5', 'NAD6', 'ATP6', 'ATP8', '16S', '12S']


def gene_of(filename):
    return filename.split('__', 1)[0]


def sample_of(filename):
    base = filename.rsplit('.', 1)[0]
    return base.split('__', 1)[1] if '__' in base else base


def read_fasta(path):
    header, seq, out = None, [], []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip('\n')
            if line.startswith('>'):
                if header is not None:
                    out.append((header, ''.join(seq)))
                header, seq = line[1:], []
            else:
                seq.append(line)
    if header is not None:
        out.append((header, ''.join(seq)))
    return out


def aggregate_ext(seqs_dir, ext, outdir):
    files = []
    if os.path.isdir(seqs_dir):
        files = sorted(f for f in os.listdir(seqs_dir) if f.endswith(ext))
    by_gene = {}
    for f in files:
        by_gene.setdefault(gene_of(f), []).append(f)
    os.makedirs(outdir, exist_ok=True)
    for stem in GENE_STEMS:
        group = by_gene.get(stem, [])
        if not group:
            continue
        records = []  # [header, seq, sample]
        for f in group:
            sample = sample_of(f)
            for header, seq in read_fasta(os.path.join(seqs_dir, f)):
                records.append([header, seq, sample])
        counts = {}
        for rec in records:
            counts[rec[0]] = counts.get(rec[0], 0) + 1
        for rec in records:
            if counts[rec[0]] > 1:
                rec[0] = '%s__%s' % (rec[0], rec[2])
        out_name = '%s.fasta' % stem if ext == '.fna' else '%s.faa' % stem
        with open(os.path.join(outdir, out_name), 'w') as fh:
            for header, seq, _ in records:
                fh.write('>%s\n%s\n' % (header, seq))


def read_occ(occ_dir):
    rows, idx = [], {}
    if not os.path.isdir(occ_dir):
        return rows, idx
    for f in sorted(os.listdir(occ_dir)):
        if not f.endswith('.occupancy.tsv'):
            continue
        with open(os.path.join(occ_dir, f)) as fh:
            header = fh.readline().rstrip('\n').split('\t')
            idx = {name: i for i, name in enumerate(header)}
            for line in fh:
                if line.strip():
                    rows.append(line.rstrip('\n').split('\t'))
    return rows, idx


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument('--seqs-dir', required=True)
    ap.add_argument('--occ-dir', required=True)
    ap.add_argument('--outdir', required=True)
    args = ap.parse_args(argv)

    aggregate_ext(args.seqs_dir, '.fna', os.path.join(args.outdir, 'nt'))
    aggregate_ext(args.seqs_dir, '.faa', os.path.join(args.outdir, 'aa'))

    rows, idx = read_occ(args.occ_dir)
    samples, data, species, flags = [], {}, {}, {}
    for parts in rows:
        if not idx or len(parts) <= idx.get('flags', 0):
            continue
        s = parts[idx['sample']]
        if s not in data:
            samples.append(s)
            data[s], flags[s] = {}, set()
            species[s] = parts[idx['species']]
        data[s][parts[idx['gene']]] = parts[idx['len_nt']]
        if parts[idx['flags']]:
            flags[s].update(parts[idx['flags']].split(','))

    os.makedirs(args.outdir, exist_ok=True)
    with open(os.path.join(args.outdir, 'occupancy.tsv'), 'w') as fh:
        fh.write('sample\tspecies\t' + '\t'.join(GENE_STEMS) + '\tflags\n')
        for s in samples:
            cells = [data[s].get(g, '0') for g in GENE_STEMS]
            fh.write('%s\t%s\t%s\t%s\n' % (s, species.get(s, ''), '\t'.join(cells),
                                           ','.join(sorted(flags[s]))))

    with open(os.path.join(args.outdir, 'mito_genes_occupancy_mqc.tsv'), 'w') as fh:
        fh.write('# id: "mito_gene_occupancy"\n')
        fh.write('# section_name: "Mito gene occupancy"\n')
        fh.write('# description: "Harvested nucleotide length per gene per sample from MITOS2 (blank = gene absent)."\n')
        fh.write('# plot_type: "table"\n')
        fh.write('Sample\t' + '\t'.join(GENE_STEMS) + '\n')
        for s in samples:
            cells = ['' if data[s].get(g, '0') == '0' else data[s].get(g, '0') for g in GENE_STEMS]
            fh.write('%s\t%s\n' % (s, '\t'.join(cells)))
    return 0


if __name__ == '__main__':
    sys.exit(main())
