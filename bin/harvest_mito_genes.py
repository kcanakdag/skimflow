#!/usr/bin/env python3
"""Harvest target mitochondrial genes from a MITOS2 annotation directory.

Reads result.fas (nucleotide) and result.faa (amino acid), normalises gene
names, keeps the 13 PCGs + 2 rRNAs, resolves duplicate copies and rRNA
fragments, and writes one single-record FASTA per gene plus a per-sample
occupancy row. Standard library only so it runs in a minimal python image.
"""
import argparse
import os
import re
import sys

# MITOS-native canonical key -> output filename stem.
FILENAME = {
    'cox1': 'COX1', 'cox2': 'COX2', 'cox3': 'COX3', 'cob': 'COB',
    'nad1': 'NAD1', 'nad2': 'NAD2', 'nad3': 'NAD3', 'nad4': 'NAD4',
    'nad4l': 'NAD4L', 'nad5': 'NAD5', 'nad6': 'NAD6',
    'atp6': 'ATP6', 'atp8': 'ATP8',
    'rrnl': '16S', 'rrns': '12S',   # rrnL -> 16S (large mt-rRNA), rrnS -> 12S (small mt-rRNA)
}
# Genes with a protein product (amino-acid output expected).
PCG = {'cox1', 'cox2', 'cox3', 'cob', 'nad1', 'nad2', 'nad3', 'nad4',
       'nad4l', 'nad5', 'nad6', 'atp6', 'atp8'}
# Fixed column / file order.
GENE_ORDER = ['cox1', 'cox2', 'cox3', 'cob', 'nad1', 'nad2', 'nad3', 'nad4',
              'nad4l', 'nad5', 'nad6', 'atp6', 'atp8', 'rrnl', 'rrns']

# Accepted synonyms (normalised: lowercase, alphanumerics only) -> canonical key.
SYNONYMS = {
    'cox1': 'cox1', 'co1': 'cox1', 'coi': 'cox1', 'coxi': 'cox1', 'mtco1': 'cox1', 'mtcox1': 'cox1',
    'cox2': 'cox2', 'co2': 'cox2', 'coii': 'cox2', 'coxii': 'cox2', 'mtco2': 'cox2', 'mtcox2': 'cox2',
    'cox3': 'cox3', 'co3': 'cox3', 'coiii': 'cox3', 'coxiii': 'cox3', 'mtco3': 'cox3', 'mtcox3': 'cox3',
    'cob': 'cob', 'cytb': 'cob', 'cyb': 'cob', 'cytochromeb': 'cob', 'mtcyb': 'cob',
    'nad1': 'nad1', 'nd1': 'nad1', 'nadh1': 'nad1',
    'nad2': 'nad2', 'nd2': 'nad2', 'nadh2': 'nad2',
    'nad3': 'nad3', 'nd3': 'nad3', 'nadh3': 'nad3',
    'nad4': 'nad4', 'nd4': 'nad4', 'nadh4': 'nad4',
    'nad4l': 'nad4l', 'nd4l': 'nad4l', 'nadh4l': 'nad4l',
    'nad5': 'nad5', 'nd5': 'nad5', 'nadh5': 'nad5',
    'nad6': 'nad6', 'nd6': 'nad6', 'nadh6': 'nad6',
    'atp6': 'atp6', 'atpase6': 'atp6', 'mtatp6': 'atp6',
    'atp8': 'atp8', 'atpase8': 'atp8', 'mtatp8': 'atp8',
    'rrnl': 'rrnl', '16s': 'rrnl', '16srrna': 'rrnl', 'lrrna': 'rrnl', 'lsurrna': 'rrnl', 'rnl': 'rrnl', 'mtrnr2': 'rrnl',
    'rrns': 'rrns', '12s': 'rrns', '12srrna': 'rrns', 'srrna': 'rrns', 'ssurrna': 'rrns', 'rns': 'rrns', 'mtrnr1': 'rrns',
}

_FRAG_RE = re.compile(r'^(?P<base>.+?)-(?P<frag>[a-z])$')              # rrnL-a, rrnL-b
_COPY_RE = re.compile(r'^(?P<base>.+?)(?:_\d+|-\d+|copy\d+)$', re.I)   # cox1_1 is what MITOS2 emits; -N/copyN are defensive


def split_tags(raw):
    """Return (base_name, frag_letter|None). A single letter AFTER a hyphen is
    an rRNA fragment tag; a trailing _N/-N/copyN is a duplicate-copy tag. Never
    strips a trailing protein letter such as the 'l' of nad4l."""
    name = raw.strip()
    m = _FRAG_RE.match(name)
    if m:
        return m.group('base'), m.group('frag')
    m = _COPY_RE.match(name)
    if m:
        return m.group('base'), None
    return name, None


def canonical_gene(raw):
    """Map a raw MITOS gene name to a canonical key, or None if not targeted."""
    base, _ = split_tags(raw)
    norm = re.sub(r'[^a-z0-9]', '', base.lower())
    return SYNONYMS.get(norm)


def parse_mitos_fasta(text):
    """Yield a dict per record: seqid,start,end,strand,raw,canon,frag,seq."""
    seqid = coords = strand = raw = None
    seq = []

    def build():
        if raw is None:
            return None
        frag = split_tags(raw)[1]
        start = end = 0
        if coords and '-' in coords:
            a, b = coords.split('-', 1)
            try:
                start, end = int(a), int(b)
            except ValueError:
                start = end = 0
        return {'seqid': seqid, 'start': start, 'end': end, 'strand': strand,
                'raw': raw, 'canon': canonical_gene(raw), 'frag': frag,
                'seq': ''.join(seq)}

    for line in text.splitlines():
        if line.startswith('>'):
            rec = build()
            if rec:
                yield rec
            fields = [f.strip() for f in line[1:].split(';')]
            seqid = fields[0] if len(fields) > 0 else ''
            coords = fields[1] if len(fields) > 1 else ''
            strand = fields[2] if len(fields) > 2 else ''
            raw = fields[3] if len(fields) > 3 else ''
            seq = []
        else:
            seq.append(line.strip())
    rec = build()
    if rec:
        yield rec


def select_features(records):
    """canon -> {seq, merged, start, end}. Merge rRNA fragments in strand-aware
    coding order; otherwise keep the longest copy."""
    by_gene = {}
    for r in records:
        if r['canon'] is not None:
            by_gene.setdefault(r['canon'], []).append(r)
    out = {}
    for canon, recs in by_gene.items():
        frags = [r for r in recs if r['frag'] is not None]
        if frags:
            # MITOS emits each fragment already in coding orientation, so coding
            # 5'->3' order is ascending genomic coordinate on the plus strand and
            # descending on the minus strand. Sort strand-aware so a minus-strand
            # rRNA is not silently scrambled (letter order is not reliable). This
            # assumes a gene is reported as fragments OR one full record, never
            # both (true for MITOS in practice).
            strand = frags[0]['strand']
            frags.sort(key=lambda r: r['start'], reverse=(strand == '-'))
            out[canon] = {'seq': ''.join(r['seq'] for r in frags), 'merged': True,
                          'start': min(r['start'] for r in frags),
                          'end': max(r['end'] for r in frags)}
        else:
            best = max(recs, key=lambda r: len(r['seq']))
            out[canon] = {'seq': best['seq'], 'merged': False,
                          'start': best['start'], 'end': best['end']}
    return out


def _overlap(a0, a1, b0, b1):
    return max(0, min(a1, b1) - max(a0, b0))


def select_aa(nt_features, aa_records):
    """For each PCG nt feature, pick the aa record (same canon) with the best
    coordinate overlap so nt and aa describe the same locus."""
    aa_by_gene = {}
    for r in aa_records:
        if r['canon'] is not None:
            aa_by_gene.setdefault(r['canon'], []).append(r)
    out = {}
    for canon, feat in nt_features.items():
        if canon not in PCG:
            continue
        cands = aa_by_gene.get(canon, [])
        if not cands:
            continue
        best = max(cands, key=lambda r: (_overlap(feat['start'], feat['end'], r['start'], r['end']),
                                         len(r['seq'])))
        out[canon] = best['seq']
    return out


def sanitize_label(s):
    """Newick-safe ASCII token: collapse whitespace to '_', drop illegal chars."""
    s = (s or '').strip()
    s = re.sub(r'\s+', '_', s)
    s = re.sub(r'[^A-Za-z0-9_.-]', '', s)
    return s or 'unknown'


def has_internal_stop(aa):
    return '*' in aa.rstrip('*')


def find_result(mitos_dir, name):
    direct = os.path.join(mitos_dir, name)
    if os.path.isfile(direct):
        return direct
    for root, _dirs, files in os.walk(mitos_dir):
        if name in files:
            return os.path.join(root, name)
    return None


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument('--mitos-dir', required=True)
    ap.add_argument('--sample-id', required=True)
    ap.add_argument('--species-id', required=True)
    ap.add_argument('--outdir', required=True)
    args = ap.parse_args(argv)

    label = sanitize_label(args.species_id) or sanitize_label(args.sample_id)
    species = label
    seqs_dir = os.path.join(args.outdir, 'seqs')
    os.makedirs(seqs_dir, exist_ok=True)

    fas = find_result(args.mitos_dir, 'result.fas')
    faa = find_result(args.mitos_dir, 'result.faa')
    nt_records = list(parse_mitos_fasta(open(fas).read())) if fas else []
    aa_records = list(parse_mitos_fasta(open(faa).read())) if faa else []

    nt_features = select_features(nt_records)
    aa_seqs = select_aa(nt_features, aa_records)

    occ = os.path.join(args.outdir, '%s.occupancy.tsv' % args.sample_id)
    with open(occ, 'w') as fh:
        fh.write('sample\tspecies\tgene\tlen_nt\tlen_aa\tflags\n')
        for canon in GENE_ORDER:
            stem = FILENAME[canon]
            feat = nt_features.get(canon)
            aa = aa_seqs.get(canon)
            len_nt = len(feat['seq']) if feat else 0
            len_aa = len(aa) if aa else 0
            flags = []
            if feat and feat['merged']:
                flags.append('merged_fragments')
            if aa and has_internal_stop(aa):
                flags.append('internal_stop')
            if feat and len_nt:
                with open(os.path.join(seqs_dir, '%s__%s.fna' % (stem, args.sample_id)), 'w') as out:
                    out.write('>%s\n%s\n' % (label, feat['seq']))
            if canon in PCG and aa:
                with open(os.path.join(seqs_dir, '%s__%s.faa' % (stem, args.sample_id)), 'w') as out:
                    out.write('>%s\n%s\n' % (label, aa))
            fh.write('%s\t%s\t%s\t%d\t%d\t%s\n' % (args.sample_id, species, stem, len_nt, len_aa, ','.join(flags)))
    return 0


if __name__ == '__main__':
    sys.exit(main())
