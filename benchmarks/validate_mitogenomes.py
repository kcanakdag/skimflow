#!/usr/bin/env python3
"""Compare skimflow mitogenomes against published GenBank references.

For every sample listed in the reference table (sample_id <TAB> accession)
this reports:

  * assembled length, number of sequences and length relative to the reference
  * blastn reference coverage (share of the reference covered by HSPs),
    assembly coverage (share of the assembly that matches the reference) and
    alignment-length-weighted identity
  * gene content of the MITOS2 annotation against the GenBank feature table,
    on the standard 37-gene set (13 PCGs, 2 rRNAs, 22 tRNAs; tRNAs compared by
    amino acid because not every record names trnL1/L2, trnS1/S2)
  * gene order: the share of reference gene adjacencies (PCGs and rRNAs,
    strand aware, circular) that the annotation reproduces

The GenBank FASTA and feature table are fetched with NCBI efetch into
--ref-dir unless they are already there. blastn must be on PATH, or pass a
command prefix, e.g. --blastn "podman run --rm -v $PWD:$PWD -w $PWD
quay.io/biocontainers/blast:2.16.0--h66d330f_5 blastn".

Stdlib only. Usage:
  benchmarks/validate_mitogenomes.py --results results --refs benchmarks/validation_refs.tsv \
      --ref-dir benchmarks/ref --out validation.tsv
"""
import argparse
import glob
import os
import re
import shlex
import subprocess
import sys
import time
import urllib.request

EFETCH = ('https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi'
          '?db=nuccore&id={acc}&rettype={rt}&retmode=text')

PCG = ['atp6', 'atp8', 'cob', 'cox1', 'cox2', 'cox3', 'nad1', 'nad2', 'nad3',
       'nad4', 'nad4l', 'nad5', 'nad6']
RRNA = ['rrnL', 'rrnS']
# 22 tRNAs as amino acid letters; L and S occur twice.
TRNA = list('ACDEFGHIKMNPQRTVWY') + ['L', 'L', 'S', 'S']

AA3 = {'Ala': 'A', 'Arg': 'R', 'Asn': 'N', 'Asp': 'D', 'Cys': 'C', 'Gln': 'Q',
       'Glu': 'E', 'Gly': 'G', 'His': 'H', 'Ile': 'I', 'Leu': 'L', 'Lys': 'K',
       'Met': 'M', 'Phe': 'F', 'Pro': 'P', 'Ser': 'S', 'Thr': 'T', 'Trp': 'W',
       'Tyr': 'Y', 'Val': 'V'}

PCG_ALIASES = {'cytb': 'cob', 'cob': 'cob', 'nd': 'nad', 'nad': 'nad',
               'co': 'cox', 'cox': 'cox', 'atp': 'atp'}


def norm_gene(name):
    """Map a GenBank or MITOS2 gene name to the standard set.

    PCGs/rRNAs -> 'cox1', 'rrnL', ...; tRNAs -> 'trnL' (amino acid only);
    anything else (OH, introns, unknown) -> None.
    """
    if not name:
        return None
    n = name.strip()
    n = re.sub(r'\(.*?\)$', '', n)        # MITOS anticodon: trnL1(tag)
    n = re.sub(r'-[a-z]$', '', n)         # MITOS split parts: cox1-a
    low = n.lower()
    m = re.fullmatch(r'trn([a-z])[12]?', low)
    if m:
        return 'trn' + m.group(1).upper()
    m = re.fullmatch(r'trna-([a-z]{3})\d?', low)
    if m and m.group(1).capitalize() in AA3:
        return 'trn' + AA3[m.group(1).capitalize()]
    if low in ('rrnl', 'rnl', 'l-rrna', '16s', '16s ribosomal rna', 'large subunit ribosomal rna'):
        return 'rrnL'
    if low in ('rrns', 'rns', 's-rrna', '12s', '12s ribosomal rna', 'small subunit ribosomal rna'):
        return 'rrnS'
    m = re.fullmatch(r'(nadh dehydrogenase|cytochrome c oxidase|atp synthase(?: f0)?) subunit (\w+)', low)
    if m:
        prefix = {'nadh dehydrogenase': 'nad', 'cytochrome c oxidase': 'cox'}.get(m.group(1), 'atp')
        num = {'i': '1', 'ii': '2', 'iii': '3'}.get(m.group(2), m.group(2))
        return prefix + num
    if low in ('cytb', 'cob', 'cyt b', 'cytochrome b'):
        return 'cob'
    m = re.fullmatch(r'(nad|nd|cox|co|atp)(\d)(l?)', low)
    if m:
        return PCG_ALIASES[m.group(1)] + m.group(2) + m.group(3)
    return None


def expected_counts():
    exp = {g: 1 for g in PCG + RRNA}
    for aa in TRNA:
        exp['trn' + aa] = exp.get('trn' + aa, 0) + 1
    return exp


def parse_feature_table(text):
    """Return [(start, end, strand, gene)] for CDS/tRNA/rRNA features.

    Names come from the gene qualifier when present, else the product.
    Coordinates are 1-based; strand '-' when start > end.
    """
    feats = []
    cur = None
    gene_at = {}
    for line in text.splitlines():
        if line.startswith('>') or not line.strip():
            continue
        cols = line.split('\t')
        if cols[0] and (len(cols) < 3 or not cols[2]) and cur is not None:
            # continuation interval of a joined feature (e.g. intron-split cox1)
            cur['e'] = int(re.sub(r'\D', '', cols[1]))
            continue
        if cols[0] and len(cols) >= 3 and cols[2]:
            s = int(re.sub(r'\D', '', cols[0]))
            e = int(re.sub(r'\D', '', cols[1]))
            cur = {'s': s, 'e': e, 'type': cols[2], 'q': {}}
            if cols[2] == 'gene':
                gene_at[(s, e)] = cur
            elif cols[2] in ('CDS', 'tRNA', 'rRNA'):
                feats.append(cur)
        elif cur is not None and len(cols) >= 5:
            cur['q'].setdefault(cols[3], cols[4])
    out = []
    gene_by_start = {s: g for (s, _), g in gene_at.items()}
    for f in feats:
        g = gene_by_start.get(f['s'], {'q': {}})['q'].get('gene')
        name = norm_gene(g) or norm_gene(f['q'].get('product'))
        if name is None:
            continue
        strand = '-' if f['s'] > f['e'] else '+'
        out.append((min(f['s'], f['e']), max(f['s'], f['e']), strand, name))
    return sorted(out)


def parse_mitos_bed(paths):
    """Return {seq_id: [(start, end, strand, gene)]} from MITOS2 result.bed files."""
    seqs = {}
    for i, p in enumerate(sorted(paths)):
        with open(p) as fh:
            for line in fh:
                cols = line.rstrip('\n').split('\t')
                if len(cols) < 6:
                    continue
                name = norm_gene(cols[3])
                if name is None:
                    continue
                key = f'{i}:{cols[0]}'
                seqs.setdefault(key, []).append((int(cols[1]), int(cols[2]), cols[5], name))
    for k in seqs:
        seqs[k].sort()
    return seqs


def mitoz_core_genes(gbk_text):
    """Distinct PCG/rRNA names in a MitoZ GenBank file (MitoZ's tRNA calls are
    not compared: they are often absent for invertebrates)."""
    names = {norm_gene(m) for m in re.findall(r'/gene="([^"]+)"', gbk_text)}
    return names & set(PCG + RRNA)


def gene_counts(features):
    counts = {}
    for f in features:
        counts[f[3]] = counts.get(f[3], 0) + 1
    return counts


def content_summary(counts):
    """(found, missing list, extra copies) against the 37-gene set."""
    exp = expected_counts()
    found = sum(min(counts.get(g, 0), n) for g, n in exp.items())
    missing = []
    for g, n in sorted(exp.items()):
        missing += [g] * max(0, n - counts.get(g, 0))
    extra = sum(max(0, c - exp.get(g, 0)) for g, c in counts.items())
    return found, missing, extra


def collapse_parts(features):
    """Merge consecutive features with the same name (split genes like cox1-a/-b)."""
    out = []
    for f in features:
        if out and out[-1][3] == f[3] and out[-1][2] == f[2]:
            continue
        out.append(f)
    return out


def adjacencies(order, circular):
    """Strand-aware adjacency set, invariant to reading the molecule backwards."""
    flip = {'+': '-', '-': '+'}
    adj = set()
    n = len(order)
    pairs = [(order[i], order[i + 1]) for i in range(n - 1)]
    if circular and n > 1:
        pairs.append((order[-1], order[0]))
    for (a, sa), (b, sb) in pairs:
        fwd = (a, sa, b, sb)
        rev = (b, flip[sb], a, flip[sa])
        adj.add(min(fwd, rev))
    return adj


def order_of(features):
    return [(f[3], f[2]) for f in collapse_parts(features) if f[3] in PCG + RRNA]


def read_fasta(path):
    seqs, name = {}, None
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if line.startswith('>'):
                name = line[1:].split()[0]
                seqs[name] = []
            elif name:
                seqs[name].append(line)
    return {k: ''.join(v) for k, v in seqs.items()}


def merge_len(intervals):
    total, cur_s, cur_e = 0, None, None
    for s, e in sorted(intervals):
        if cur_e is None or s > cur_e + 1:
            if cur_e is not None:
                total += cur_e - cur_s + 1
            cur_s, cur_e = s, e
        else:
            cur_e = max(cur_e, e)
    if cur_e is not None:
        total += cur_e - cur_s + 1
    return total


def blast_stats(hsps, qlen, slen):
    """hsps: [(pident, length, qstart, qend, sstart, send, bitscore, qseqid)].

    Identity is weighted over a greedy, bitscore-ordered set of HSPs that do
    not overlap on the reference by more than half, so repeats in the
    assembly do not count twice.
    """
    chosen, taken = [], []
    for h in sorted(hsps, key=lambda h: -h[6]):
        s, e = sorted((h[4], h[5]))
        ov = sum(max(0, min(e, te) - max(s, ts) + 1) for ts, te in taken)
        if ov > 0.5 * (e - s + 1):
            continue
        chosen.append(h)
        taken.append((s, e))
    if not chosen:
        return 0.0, 0.0, None
    ref_cov = merge_len([sorted((h[4], h[5])) for h in hsps]) / slen
    by_seq = {}
    for h in hsps:
        by_seq.setdefault(h[7], []).append(sorted((h[2], h[3])))
    asm_cov = sum(merge_len(iv) for iv in by_seq.values()) / qlen
    ident = sum(h[0] * h[1] for h in chosen) / sum(h[1] for h in chosen)
    return ref_cov, asm_cov, ident


def run_blast(blastn, query, subject):
    fmt = '6 qseqid pident length qstart qend sstart send bitscore evalue'
    cmd = shlex.split(blastn) + ['-query', query, '-subject', subject,
                                 '-outfmt', fmt, '-evalue', '1e-10']
    out = subprocess.run(cmd, check=True, capture_output=True, text=True).stdout
    hsps = []
    for line in out.splitlines():
        c = line.split('\t')
        if int(c[2]) < 100:
            continue
        hsps.append((float(c[1]), int(c[2]), int(c[3]), int(c[4]),
                     int(c[5]), int(c[6]), float(c[7]), c[0]))
    return hsps


def fetch(acc, ref_dir):
    os.makedirs(ref_dir, exist_ok=True)
    for rt, ext in (('fasta', 'fasta'), ('ft', 'ft')):
        dest = os.path.join(ref_dir, f'{acc}.{ext}')
        if os.path.exists(dest) and os.path.getsize(dest) > 0:
            continue
        with urllib.request.urlopen(EFETCH.format(acc=acc, rt=rt), timeout=60) as r:
            data = r.read().decode()
        with open(dest, 'w') as fh:
            fh.write(data)
        time.sleep(0.4)  # NCBI asks for <= 3 requests/s without an API key


def validate(sample, acc, results, ref_dir, blastn):
    row = {'sample_id': sample, 'accession': acc}
    ref = read_fasta(os.path.join(ref_dir, f'{acc}.fasta'))
    slen = sum(len(s) for s in ref.values())
    with open(os.path.join(ref_dir, f'{acc}.ft')) as fh:
        ref_feats = parse_feature_table(fh.read())
    row['ref_bp'] = slen

    asm_path = os.path.join(results, 'mitogenome', sample, f'{sample}.mito.fasta')
    if not os.path.exists(asm_path):
        row['status'] = 'no_assembly'
        return row
    asm = read_fasta(asm_path)
    qlen = sum(len(s) for s in asm.values())
    row.update(status='ok', asm_bp=qlen, asm_seqs=len(asm),
               length_ratio=round(qlen / slen, 3))

    hsps = run_blast(blastn, asm_path, os.path.join(ref_dir, f'{acc}.fasta'))
    ref_cov, asm_cov, ident = blast_stats(hsps, qlen, slen)
    row.update(ref_cov_pct=round(100 * ref_cov, 1), asm_cov_pct=round(100 * asm_cov, 1),
               identity_pct=round(ident, 2) if ident is not None else '')

    ref_found, _, _ = content_summary(gene_counts(ref_feats))
    row['ref_genes'] = ref_found
    gbk = os.path.join(results, 'mitogenome_annotation', 'mitoz', sample, f'{sample}.mitoz.gbk')
    if os.path.exists(gbk):
        with open(gbk) as fh:
            core = mitoz_core_genes(fh.read())
        row['mitoz_core'] = f'{len(core)}/{len(PCG + RRNA)}'
    beds = glob.glob(os.path.join(results, 'mitogenome_annotation', 'mitos2', sample,
                                  '**', 'result.bed'), recursive=True)
    if not beds:
        row['mitos2'] = 'missing'
        return row
    ann = parse_mitos_bed(beds)
    all_feats = [f for fs in ann.values() for f in fs]
    found, missing, extra = content_summary(
        gene_counts([f for fs in ann.values() for f in collapse_parts(fs)]))
    row.update(mitos2_genes=found, missing=','.join(missing) or '-', extra_copies=extra)

    ref_adj = adjacencies(order_of(ref_feats), circular=True)
    circ = len(ann) == 1 and 'circular' in next(iter(ann)).lower()
    asm_adj = set()
    for fs in ann.values():
        asm_adj |= adjacencies(order_of(fs), circular=circ)
    row['order_adj'] = f'{len(ref_adj & asm_adj)}/{len(ref_adj)}'
    row['annotated_features'] = len(all_feats)
    return row


COLUMNS = ['sample_id', 'accession', 'status', 'ref_bp', 'asm_bp', 'asm_seqs',
           'length_ratio', 'ref_cov_pct', 'asm_cov_pct', 'identity_pct',
           'ref_genes', 'mitos2_genes', 'missing', 'extra_copies', 'order_adj',
           'mitoz_core']


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('--results', required=True, help='pipeline outdir')
    ap.add_argument('--refs', required=True, help='TSV: sample_id, accession')
    ap.add_argument('--ref-dir', required=True, help='cache for GenBank FASTA + feature tables')
    ap.add_argument('--blastn', default='blastn', help='blastn command (prefix)')
    ap.add_argument('--out', default='-', help='output TSV (default stdout)')
    a = ap.parse_args(argv)

    pairs = []
    with open(a.refs) as fh:
        for line in fh:
            if line.startswith('#') or not line.strip():
                continue
            s, acc = line.split()[:2]
            if s == 'sample_id':
                continue
            pairs.append((s, acc))

    rows = []
    for s, acc in pairs:
        fetch(acc, a.ref_dir)
        rows.append(validate(s, acc, a.results, a.ref_dir, a.blastn))

    out = sys.stdout if a.out == '-' else open(a.out, 'w')
    out.write('\t'.join(COLUMNS) + '\n')
    for r in rows:
        out.write('\t'.join(str(r.get(c, '')) for c in COLUMNS) + '\n')
    if out is not sys.stdout:
        out.close()


if __name__ == '__main__':
    main()
