import os
import subprocess
import sys
import tempfile
import unittest

BIN = os.path.join(os.path.dirname(__file__), '..', 'bin')
SCRIPT = os.path.join(BIN, 'concat_alignments.py')


def _write(path, recs):
    with open(path, 'w') as fh:
        for name, seq in recs:
            fh.write('>%s\n%s\n' % (name, seq))


def _read_fasta(path):
    out, name = {}, None
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if line.startswith('>'):
                name = line[1:]
                out[name] = ''
            elif name:
                out[name] += line
    return out


class TestConcat(unittest.TestCase):
    def _run(self, min_taxa):
        d = tempfile.mkdtemp()
        aln = os.path.join(d, 'aln')
        out = os.path.join(d, 'out')
        os.makedirs(aln)
        # 16S given first on disk to check canonical gene order (COX1 before 16S).
        _write(os.path.join(aln, '16S.trimmed.fasta'),
               [('Taxon_A', 'GGGG'), ('Taxon_B', 'GG-G')])
        _write(os.path.join(aln, 'COX1.trimmed.fasta'),
               [('Taxon_A', 'ATGAAA'), ('Taxon_B', 'ATGA-A'), ('Taxon_C', 'ATGCCC')])
        res = subprocess.run(
            [sys.executable, SCRIPT, '--aln-dir', aln, '--outdir', out,
             '--min-taxa', str(min_taxa)],
            check=True, capture_output=True, text=True)
        return out, res.stderr

    def test_missing_gene_padded_with_gaps(self):
        out, _ = self._run(2)
        sm = _read_fasta(os.path.join(out, 'supermatrix.fasta'))
        self.assertEqual(sorted(sm), ['Taxon_A', 'Taxon_B', 'Taxon_C'])
        self.assertEqual(sm['Taxon_A'], 'ATGAAA' + 'GGGG')
        self.assertEqual(sm['Taxon_C'], 'ATGCCC' + '----')
        self.assertEqual({len(s) for s in sm.values()}, {10})

    def test_partition_coordinates(self):
        out, _ = self._run(2)
        with open(os.path.join(out, 'partitions.nex')) as fh:
            text = fh.read()
        self.assertIn('charset COX1 = 1-6;', text)
        self.assertIn('charset 16S = 7-10;', text)
        self.assertLess(text.index('COX1'), text.index('16S'))
        with open(os.path.join(out, 'concat_stats.tsv')) as fh:
            lines = fh.read().splitlines()
        self.assertEqual(lines[1], '3\t2\t10\tCOX1,16S')

    def test_gene_below_min_taxa_dropped(self):
        out, err = self._run(3)
        sm = _read_fasta(os.path.join(out, 'supermatrix.fasta'))
        self.assertEqual(sm['Taxon_A'], 'ATGAAA')
        with open(os.path.join(out, 'partitions.nex')) as fh:
            text = fh.read()
        self.assertNotIn('16S', text)
        self.assertIn('16S: 2 taxa < --min-taxa 3; dropped', err)


if __name__ == '__main__':
    unittest.main()
