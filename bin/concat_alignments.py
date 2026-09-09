#!/usr/bin/env python3
"""Concatenate per-gene trimmed alignments into one supermatrix + a partition
file, for a single sequence type (nt or aa).

Each input alignment is a FASTA whose headers are species/sample labels. A
taxon missing a gene is padded with gaps for that gene's alignment length, so
every taxon spans the full supermatrix (standard supermatrix / total-evidence
concatenation). Genes present in fewer than --min-taxa sequences are dropped
(too sparse to inform the tree). Standard library only so it runs in a minimal
python image.

Outputs (into --outdir):
  supermatrix.fasta   concatenated alignment (one record per taxon)
  partitions.nex      NEXUS charset block, one charset per gene (IQ-TREE -p)
  concat_stats.tsv    n_taxa, n_genes, alignment_length, genes_used
"""
import argparse
import os
import sys

# Canonical gene column / partition order (matches bin/aggregate_mito_genes.py).
GENE_ORDER = ['COX1', 'COX2', 'COX3', 'COB', 'NAD1', 'NAD2', 'NAD3', 'NAD4',
              'NAD4L', 'NAD5', 'NAD6', 'ATP6', 'ATP8', '16S', '12S']


def read_fasta(path):
    header, seq, out = None, [], []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip('\n')
            if line.startswith('>'):
                if header is not None:
                    out.append((header, ''.join(seq)))
                header, seq = line[1:].strip(), []
            else:
                seq.append(line.strip())
    if header is not None:
        out.append((header, ''.join(seq)))
    return out


def gene_of(filename):
    """Gene stem is the filename up to the first dot: COX1.trimmed.fasta -> COX1."""
    return os.path.basename(filename).split('.', 1)[0]


def gene_sort_key(stem):
    return GENE_ORDER.index(stem) if stem in GENE_ORDER else len(GENE_ORDER)


def load_genes(aln_dir):
    """gene_stem -> list[(taxon, seq)], only for non-empty alignment files."""
    genes = {}
    if not os.path.isdir(aln_dir):
        return genes
    for f in sorted(os.listdir(aln_dir)):
        path = os.path.join(aln_dir, f)
        if not os.path.isfile(path):
            continue
        recs = [(h, s) for h, s in read_fasta(path) if s]
        if recs:
            genes[gene_of(f)] = recs
    return genes


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument('--aln-dir', required=True,
                    help='directory of per-gene trimmed alignment FASTAs')
    ap.add_argument('--outdir', required=True)
    ap.add_argument('--min-taxa', type=int, default=4,
                    help='drop a gene present in fewer than this many taxa')
    args = ap.parse_args(argv)

    os.makedirs(args.outdir, exist_ok=True)
    genes = load_genes(args.aln_dir)

    # Keep genes in canonical order, each with a single (equal) alignment width.
    used = []       # [(stem, length, {taxon: seq})]
    for stem in sorted(genes, key=gene_sort_key):
        recs = genes[stem]
        if len({len(s) for _, s in recs}) != 1:
            sys.stderr.write(
                '[concat] %s: sequences are not equal length (not aligned); '
                'skipping\n' % stem)
            continue
        if len(recs) < args.min_taxa:
            continue
        length = len(recs[0][1])
        # If a taxon appears twice for one gene (should not happen post-harvest),
        # keep the first occurrence deterministically.
        by_taxon = {}
        for taxon, seq in recs:
            by_taxon.setdefault(taxon, seq)
        used.append((stem, length, by_taxon))

    taxa = sorted({t for _, _, bt in used for t in bt})

    # Build the concatenated rows; missing gene -> gap run of that gene's width.
    matrix = {t: [] for t in taxa}
    for _stem, length, by_taxon in used:
        gap = '-' * length
        for t in taxa:
            matrix[t].append(by_taxon.get(t, gap))

    supermatrix = os.path.join(args.outdir, 'supermatrix.fasta')
    with open(supermatrix, 'w') as fh:
        for t in taxa:
            fh.write('>%s\n%s\n' % (t, ''.join(matrix[t])))

    # NEXUS partition file: one charset per gene. No datatype/model token, so
    # IQ-TREE auto-detects nt vs aa from the alignment and ModelFinder (-m MFP)
    # picks a model per partition.
    partitions = os.path.join(args.outdir, 'partitions.nex')
    with open(partitions, 'w') as fh:
        fh.write('#nexus\nbegin sets;\n')
        pos = 1
        for stem, length, _bt in used:
            end = pos + length - 1
            fh.write('  charset %s = %d-%d;\n' % (stem, pos, end))
            pos = end + 1
        fh.write('end;\n')

    total_len = sum(length for _s, length, _bt in used)
    stats = os.path.join(args.outdir, 'concat_stats.tsv')
    with open(stats, 'w') as fh:
        fh.write('n_taxa\tn_genes\talignment_length\tgenes_used\n')
        fh.write('%d\t%d\t%d\t%s\n' % (
            len(taxa), len(used), total_len,
            ','.join(stem for stem, _l, _bt in used)))

    sys.stderr.write('[concat] %d taxa x %d genes = %d columns\n'
                     % (len(taxa), len(used), total_len))
    return 0


if __name__ == '__main__':
    sys.exit(main())
