#!/usr/bin/env python3
"""Pass/fail check for the skimflow benchmark.

Reads the table written by validate_mitogenomes.py and compares every sample
against the thresholds in expected.tsv (one row per sample). Optionally checks
that the species trees contain the expected splits (expected_splits.tsv: tree
name, then a comma-separated taxon set). Trees are unrooted, so a split
matches either side of a branch.

Prints one line per check and exits 0 when every check passes, 1 otherwise.

Stdlib only. Usage:
  benchmarks/check_benchmark.py --validation results_benchmark/validation.tsv \\
      --expected benchmarks/expected.tsv \\
      --phylogeny results_benchmark/phylogeny --splits benchmarks/expected_splits.tsv
"""
import argparse
import csv
import os
import re
import sys

# (column in validation.tsv, threshold column in expected.tsv, comparison)
CHECKS = [
    ('ref_cov_pct', 'min_ref_cov_pct', '>='),
    ('identity_pct', 'min_identity_pct', '>='),
    ('mitos2_genes', 'min_mitos2_genes', '>='),
    ('asm_seqs', 'max_asm_seqs', '<='),
]


def read_tsv(path):
    with open(path, newline='') as fh:
        return list(csv.DictReader((l for l in fh if not l.startswith('#')), delimiter='\t'))


def _num(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return None


def check_samples(validation_path, expected_path):
    """Return a list of (sample, check name, ok, message)."""
    got = {r['sample_id']: r for r in read_tsv(validation_path)}
    out = []
    for exp in read_tsv(expected_path):
        s = exp['sample_id']
        r = got.get(s)
        if r is None or r.get('status') != 'ok':
            status = 'not in validation table' if r is None else r.get('status')
            out.append((s, 'status', False, f'no usable mitogenome ({status})'))
            continue
        for col, limit_col, op in CHECKS:
            limit = _num(exp.get(limit_col))
            if limit is None:
                continue
            value = _num(r.get(col))
            ok = value is not None and (value >= limit if op == '>=' else value <= limit)
            out.append((s, col, ok, f'{r.get(col) or "NA"} (need {op} {exp[limit_col]})'))
        lo, hi = _num(exp.get('min_length_ratio')), _num(exp.get('max_length_ratio'))
        if lo is not None or hi is not None:
            v = _num(r.get('length_ratio'))
            ok = v is not None and (lo is None or v >= lo) and (hi is None or v <= hi)
            out.append((s, 'length_ratio', ok,
                        f'{r.get("length_ratio") or "NA"} (need {lo} to {hi})'))
    return out


def parse_newick(text):
    """Return (leaf names, list of clades as frozensets of leaf names)."""
    text = text.strip().rstrip(';')
    leaves, clades, stack = set(), [], []
    i, n = 0, len(text)
    while i < n:
        c = text[i]
        if c == '(':
            stack.append(set())
            i += 1
        elif c == ',':
            i += 1
        elif c == ')':
            members = stack.pop()
            clades.append(frozenset(members))
            if stack:
                stack[-1] |= members
            i += 1
            # skip the support label and branch length of this clade
            m = re.match(r"[^,():;]*(:[^,();]*)?", text[i:])
            i += m.end()
        else:
            if c == "'":
                j = text.index("'", i + 1)
                name = text[i + 1:j]
                i = j + 1
            else:
                m = re.match(r"[^,():;]+", text[i:])
                if m is None:
                    raise ValueError(f'bad Newick near position {i}: {text[i:i + 20]!r}')
                name = m.group(0).strip()
                i += m.end()
            m = re.match(r"(:[^,();]*)?", text[i:])
            i += m.end()
            leaves.add(name)
            if stack:
                stack[-1].add(name)
    return leaves, clades


def has_split(newick, taxa):
    leaves, clades = parse_newick(newick)
    taxa = frozenset(taxa)
    other = frozenset(leaves) - taxa
    return any(c == taxa or c == other for c in clades)


def check_splits(phylo_dir, splits_path):
    out = []
    for row in read_tsv(splits_path):
        tree, taxa = row['tree'], [t.strip() for t in row['taxa'].split(',') if t.strip()]
        path = os.path.join(phylo_dir, tree, f'{tree}_species.treefile')
        label = f'{tree} tree'
        if not os.path.exists(path):
            out.append((label, '{' + ','.join(taxa) + '}', False, f'missing {path}'))
            continue
        with open(path) as fh:
            newick = fh.read()
        leaves, _ = parse_newick(newick)
        absent = [t for t in taxa if t not in leaves]
        if absent:
            out.append((label, '{' + ','.join(taxa) + '}', False,
                        f'taxa not in tree: {",".join(absent)}'))
            continue
        ok = has_split(newick, taxa)
        out.append((label, '{' + ','.join(taxa) + '}', ok,
                    'split present' if ok else 'split not found'))
    return out


def main(argv=None):
    ap = argparse.ArgumentParser(description='Pass/fail check for the skimflow benchmark.')
    ap.add_argument('--validation', required=True, help='validate_mitogenomes.py output')
    ap.add_argument('--expected', required=True, help='per-sample thresholds TSV')
    ap.add_argument('--phylogeny', help='pipeline phylogeny dir (results/phylogeny)')
    ap.add_argument('--splits', help='expected tree splits TSV')
    args = ap.parse_args(argv)

    results = check_samples(args.validation, args.expected)
    if args.phylogeny and args.splits:
        results += check_splits(args.phylogeny, args.splits)

    width = max(len(s) for s, *_ in results) if results else 10
    for subject, name, ok, msg in results:
        print(f'{"PASS" if ok else "FAIL"}  {subject:<{width}}  {name:<14} {msg}')
    failed = sum(1 for *_, ok, _ in results if not ok)
    print(f'\n{len(results) - failed}/{len(results)} checks passed'
          + ('' if failed == 0 else f', {failed} failed'))
    return 0 if failed == 0 else 1


if __name__ == '__main__':
    sys.exit(main())
