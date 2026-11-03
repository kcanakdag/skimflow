#!/usr/bin/env python3
"""Pick a GetOrganelle seed from a sample's own assembly ("self-seeding").

Input is a tblastn table of mitochondrial protein-coding gene proteins
(query ids "<gene>|<species>|<accession>") searched against the sample's
contigs with the mitochondrial genetic code, plus the contig FASTA.

A contig is a candidate when it has at least one hit that covers enough of
its query protein (--min-qcov) and has no stop codon in the translated
contig sequence (stops suggest a NUMT, a decayed nuclear copy, or the wrong
genetic code). Mitochondrial contigs have far higher read depth than nuclear
ones, and a contaminant's mitogenome usually has lower depth than the
sample's own, so the deepest candidate anchors the selection: every
candidate with at least --depth-frac of its depth is kept, deepest first,
until --max-seed-bp, as long as it adds a gene the seed does not have yet
(a contig that only repeats genes is a second copy: another organism's
mitogenome or a nuclear copy). Depth is MEGAHIT's multi= value in the contig header;
without it, candidates are ranked by bitscore and all are kept.

Writes the seed FASTA (whole contigs; only when something was found) and a
one-row MultiQC table. Standard library only.
"""
import argparse
import re
import statistics
import sys

HIT_FIELDS = ['qseqid', 'sseqid', 'pident', 'length', 'qstart', 'qend',
              'sstart', 'send', 'evalue', 'bitscore', 'qlen', 'sseq']


def read_contigs(path):
    """{name: (sequence, depth or None)}; depth from MEGAHIT 'multi='."""
    contigs, name, depth, chunks = {}, None, None, []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip('\n')
            if line.startswith('>'):
                if name is not None:
                    contigs[name] = (''.join(chunks), depth)
                head = line[1:].split()
                name = head[0] if head else ''
                m = re.search(r'\bmulti=([0-9.eE+-]+)', line)
                depth = float(m.group(1)) if m else None
                chunks = []
            elif name is not None:
                chunks.append(line.strip())
    if name is not None:
        contigs[name] = (''.join(chunks), depth)
    return contigs


def read_hits(path):
    hits = []
    with open(path) as fh:
        for line in fh:
            if not line.strip() or line.startswith('#'):
                continue
            c = line.rstrip('\n').split('\t')
            if len(c) < len(HIT_FIELDS):
                continue
            h = dict(zip(HIT_FIELDS, c))
            for k in ('qstart', 'qend', 'sstart', 'send', 'qlen', 'length'):
                h[k] = int(h[k])
            h['bitscore'] = float(h['bitscore'])
            h['gene'] = h['qseqid'].split('|')[0]
            hits.append(h)
    return hits


def candidates(hits, min_qcov):
    """{contig: {'genes': set, 'bitscore': best}} and rejection counts."""
    cand, rejected = {}, {'stop_codon': 0, 'partial': 0}
    for h in hits:
        if '*' in h['sseq']:
            rejected['stop_codon'] += 1
            continue
        if (h['qend'] - h['qstart'] + 1) / h['qlen'] < min_qcov:
            rejected['partial'] += 1
            continue
        c = cand.setdefault(h['sseqid'], {'genes': set(), 'bitscore': 0.0})
        c['genes'].add(h['gene'])
        c['bitscore'] = max(c['bitscore'], h['bitscore'])
    return cand, rejected


def select(cand, contigs, depth_frac, max_seed_bp):
    """Return (kept contig names, note)."""
    known = {n: contigs[n][1] for n in cand if n in contigs}
    if not known:
        return [], 'no candidate contig'
    have_depth = all(d is not None for d in known.values())
    if have_depth:
        order = sorted(known, key=lambda n: (-known[n], -cand[n]['bitscore']))
        floor = depth_frac * known[order[0]]
        order = [n for n in order if known[n] >= floor]
    else:
        order = sorted(known, key=lambda n: -cand[n]['bitscore'])
    kept, total, covered, redundant = [], 0, set(), []
    for n in order:
        ln = len(contigs[n][0])
        if kept and total + ln > max_seed_bp:
            break
        if kept and cand[n]['genes'] <= covered:
            # Only genes already in the seed: a second copy (another
            # organism's mitogenome or a nuclear copy), not a missing piece.
            redundant.append(n)
            continue
        kept.append(n)
        total += ln
        covered |= cand[n]['genes']
    note = []
    if redundant:
        note.append(f'{len(redundant)} contig(s) with only duplicate genes dropped '
                    f'(second mitogenome or NUMT?): {",".join(redundant)}')
    if not have_depth:
        note.append('no depth in headers; ranked by score')
    best_score = max(cand, key=lambda n: cand[n]['bitscore'])
    if best_score in known and best_score not in kept:
        note.append(f'best-scoring contig {best_score} not kept (low depth: contamination or NUMT?)')
    return kept, '; '.join(note) or '-'


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('--hits', required=True)
    ap.add_argument('--contigs', required=True)
    ap.add_argument('--sample', required=True)
    ap.add_argument('--out-fasta', required=True)
    ap.add_argument('--out-tsv', required=True)
    ap.add_argument('--min-qcov', type=float, default=0.3)
    ap.add_argument('--depth-frac', type=float, default=0.2)
    ap.add_argument('--max-seed-bp', type=int, default=30000)
    a = ap.parse_args(argv)

    contigs = read_contigs(a.contigs)
    cand, rejected = candidates(read_hits(a.hits), a.min_qcov)
    kept, note = select(cand, contigs, a.depth_frac, a.max_seed_bp)

    genes = sorted({g for n in kept for g in cand[n]['genes']})
    depths = [d for _, d in contigs.values() if d is not None]
    median = statistics.median(depths) if depths else None
    top = contigs[kept[0]][1] if kept else None
    seed_bp = sum(len(contigs[n][0]) for n in kept)

    if kept:
        with open(a.out_fasta, 'w') as fh:
            for n in kept:
                seq = contigs[n][0]
                fh.write(f'>{n}\n')
                for i in range(0, len(seq), 80):
                    fh.write(seq[i:i + 80] + '\n')

    fmt = lambda x: '' if x is None else f'{x:.1f}'
    with open(a.out_tsv, 'w') as fh:
        fh.write('# id: "mito_self_seed"\n'
                 '# section_name: "Mitogenome seed search"\n'
                 '# description: "Seed for GetOrganelle taken from the sample\'s own '
                 'MEGAHIT contigs (tblastn with mitochondrial proteins; deepest '
                 'clean hits kept)."\n'
                 '# plot_type: "table"\n')
        fh.write('Sample\tStatus\tSeed_contigs\tSeed_bp\tGenes_hit\tTop_depth\t'
                 'Median_contig_depth\tCandidates\tRejected_stop\tRejected_partial\tNote\n')
        fh.write('\t'.join([
            a.sample, 'found' if kept else 'none', str(len(kept)), str(seed_bp),
            f'{len(genes)}/13', fmt(top), fmt(median), str(len(cand)),
            str(rejected['stop_codon']), str(rejected['partial']), note,
        ]) + '\n')


if __name__ == '__main__':
    sys.exit(main())
