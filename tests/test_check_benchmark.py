import os
import sys
import tempfile
import unittest

BENCH = os.path.join(os.path.dirname(__file__), '..', 'benchmarks')
sys.path.insert(0, BENCH)
import check_benchmark as cb  # noqa: E402

HEADER = ('sample_id\taccession\tstatus\tref_bp\tasm_bp\tasm_seqs\tlength_ratio\t'
          'ref_cov_pct\tasm_cov_pct\tidentity_pct\tref_genes\tmitos2_genes\tmissing\t'
          'extra_copies\torder_adj\tmitoz_core\n')
EXPECTED = ('sample_id\tmin_ref_cov_pct\tmin_identity_pct\tmin_mitos2_genes\t'
            'max_asm_seqs\tmin_length_ratio\tmax_length_ratio\n'
            'A\t95\t99.5\t35\t3\t0.9\t1.25\n'
            'B\t95\t99.5\t35\t3\t0.9\t1.25\n')


def row(sample, status='ok', cov='98.0', ident='99.9', genes='37', seqs='1', ratio='1.0'):
    return (f'{sample}\tX\t{status}\t100\t100\t{seqs}\t{ratio}\t{cov}\t100\t{ident}\t'
            f'37\t{genes}\t-\t0\t15/15\t15/15\n')


class TestNewick(unittest.TestCase):
    def test_leaves_and_clades(self):
        leaves, clades = cb.parse_newick('(A:0.1,((B:1,C:1)90/80:0.2,D:1)100/100:0.1,E:0.2);')
        self.assertEqual(leaves, {'A', 'B', 'C', 'D', 'E'})
        self.assertIn(frozenset({'B', 'C'}), clades)
        self.assertIn(frozenset({'B', 'C', 'D'}), clades)

    def test_unrooted_split_matches_complement(self):
        tree = '(A:0.1,((B:1,C:1)90/80:0.2,D:1)100/100:0.1,E:0.2);'
        # {A,E} is not a clade in this rooting, but it is the complement of {B,C,D}.
        self.assertTrue(cb.has_split(tree, {'A', 'E'}))
        self.assertTrue(cb.has_split(tree, {'B', 'C'}))
        self.assertFalse(cb.has_split(tree, {'A', 'B'}))

    def test_quoted_and_unlabelled(self):
        self.assertTrue(cb.has_split("(('x y',b),(c,d));", {"x y", 'b'}))


class TestChecks(unittest.TestCase):
    def _files(self, validation_rows):
        d = tempfile.mkdtemp()
        v = os.path.join(d, 'validation.tsv')
        e = os.path.join(d, 'expected.tsv')
        with open(v, 'w') as fh:
            fh.write(HEADER + ''.join(validation_rows))
        with open(e, 'w') as fh:
            fh.write(EXPECTED)
        return v, e

    def test_all_pass(self):
        v, e = self._files([row('A'), row('B')])
        results = cb.check_samples(v, e)
        self.assertTrue(all(ok for _, _, ok, _ in results))

    def test_low_coverage_fails(self):
        v, e = self._files([row('A', cov='80.0'), row('B')])
        failed = [(s, name) for s, name, ok, _ in cb.check_samples(v, e) if not ok]
        self.assertEqual(failed, [('A', 'ref_cov_pct')])

    def test_missing_sample_fails(self):
        v, e = self._files([row('A')])
        failed = [(s, name) for s, name, ok, _ in cb.check_samples(v, e) if not ok]
        self.assertEqual(failed, [('B', 'status')])

    def test_no_assembly_fails(self):
        v, e = self._files([row('A'), 'B\tX\tno_assembly\t100\n'])
        failed = [(s, name) for s, name, ok, _ in cb.check_samples(v, e) if not ok]
        self.assertEqual(failed, [('B', 'status')])

    def test_length_ratio_bounds(self):
        v, e = self._files([row('A', ratio='1.4'), row('B', ratio='0.5')])
        failed = sorted((s, name) for s, name, ok, _ in cb.check_samples(v, e) if not ok)
        self.assertEqual(failed, [('A', 'length_ratio'), ('B', 'length_ratio')])


if __name__ == '__main__':
    unittest.main()
