import os
import subprocess
import sys
import tempfile
import unittest

BIN = os.path.join(os.path.dirname(__file__), '..', 'bin')
sys.path.insert(0, BIN)
import aggregate_mito_genes as a  # noqa: E402


def write(path, text):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, 'w') as fh:
        fh.write(text)


class TestHelpers(unittest.TestCase):
    def test_gene_and_sample_of(self):
        self.assertEqual(a.gene_of('COX1__OkadaiA.fna'), 'COX1')
        self.assertEqual(a.sample_of('COX1__OkadaiA.fna'), 'OkadaiA')
        self.assertEqual(a.gene_of('16S__TuberB.faa'), '16S')


class TestAggregate(unittest.TestCase):
    def _build(self, root):
        seqs = os.path.join(root, 'seqs')
        occ = os.path.join(root, 'occ')
        # two samples, same species name -> header collision must disambiguate
        write(os.path.join(seqs, 'COX1__OkadaiA.fna'), '>Syllis_okadai\nAAAA\n')
        write(os.path.join(seqs, 'COX1__OkadaiB.fna'), '>Syllis_okadai\nCCCC\n')
        write(os.path.join(seqs, 'COX1__OkadaiA.faa'), '>Syllis_okadai\nMK\n')
        write(os.path.join(seqs, 'COX1__OkadaiB.faa'), '>Syllis_okadai\nML\n')
        write(os.path.join(occ, 'OkadaiA.occupancy.tsv'),
              'sample\tspecies\tgene\tlen_nt\tlen_aa\tflags\n'
              'OkadaiA\tSyllis_okadai\tCOX1\t4\t2\t\n'
              'OkadaiA\tSyllis_okadai\t16S\t0\t0\t\n')
        write(os.path.join(occ, 'OkadaiB.occupancy.tsv'),
              'sample\tspecies\tgene\tlen_nt\tlen_aa\tflags\n'
              'OkadaiB\tSyllis_okadai\tCOX1\t4\t2\t\n'
              'OkadaiB\tSyllis_okadai\t16S\t0\t0\t\n')
        return seqs, occ

    def test_cli_outputs(self):
        with tempfile.TemporaryDirectory() as root:
            seqs, occ = self._build(root)
            out = os.path.join(root, 'out')
            subprocess.run(
                [sys.executable, os.path.join(BIN, 'aggregate_mito_genes.py'),
                 '--seqs-dir', seqs, '--occ-dir', occ, '--outdir', out],
                check=True)
            with open(os.path.join(out, 'nt', 'COX1.fasta')) as fh:
                cox1 = fh.read()
            # both records present, headers disambiguated by sample on collision
            self.assertEqual(cox1.count('>'), 2)
            self.assertIn('>Syllis_okadai__OkadaiA', cox1)
            self.assertIn('>Syllis_okadai__OkadaiB', cox1)
            self.assertTrue(os.path.exists(os.path.join(out, 'aa', 'COX1.faa')))
            # 16S absent in every sample -> no file emitted
            self.assertFalse(os.path.exists(os.path.join(out, 'nt', '16S.fasta')))
            with open(os.path.join(out, 'occupancy.tsv')) as fh:
                matrix = fh.read()
            self.assertIn('COX1', matrix.splitlines()[0])
            with open(os.path.join(out, 'mito_genes_occupancy_mqc.tsv')) as fh:
                mqc = fh.read()
            self.assertIn('plot_type: "table"', mqc)
            self.assertIn('# section_name:', mqc)

    def test_empty_inputs_do_not_crash(self):
        with tempfile.TemporaryDirectory() as root:
            out = os.path.join(root, 'out')
            subprocess.run(
                [sys.executable, os.path.join(BIN, 'aggregate_mito_genes.py'),
                 '--seqs-dir', os.path.join(root, 'none'),
                 '--occ-dir', os.path.join(root, 'none'),
                 '--outdir', out],
                check=True)
            self.assertTrue(os.path.exists(os.path.join(out, 'occupancy.tsv')))
            self.assertTrue(os.path.exists(os.path.join(out, 'mito_genes_occupancy_mqc.tsv')))

    def test_asymmetric_nt_aa_collision(self):
        """Regression: sample A has nt+aa, sample B has nt only.
        Both share species header >Syllis_okadai for COX1.
        Both nt and aa outputs must suffix sample A's header consistently."""
        with tempfile.TemporaryDirectory() as root:
            seqs = os.path.join(root, 'seqs')
            occ = os.path.join(root, 'occ')
            # sample A: nt and aa present
            write(os.path.join(seqs, 'COX1__A.fna'), '>Syllis_okadai\nAAAA\n')
            write(os.path.join(seqs, 'COX1__A.faa'), '>Syllis_okadai\nMK\n')
            # sample B: nt only, no aa
            write(os.path.join(seqs, 'COX1__B.fna'), '>Syllis_okadai\nCCCC\n')
            write(os.path.join(occ, 'dummy.occupancy.tsv'),
                  'sample\tspecies\tgene\tlen_nt\tlen_aa\tflags\n')
            out = os.path.join(root, 'out')
            subprocess.run(
                [sys.executable, os.path.join(BIN, 'aggregate_mito_genes.py'),
                 '--seqs-dir', seqs, '--occ-dir', occ, '--outdir', out],
                check=True)
            with open(os.path.join(out, 'nt', 'COX1.fasta')) as fh:
                nt = fh.read()
            with open(os.path.join(out, 'aa', 'COX1.faa')) as fh:
                aa = fh.read()
            # nt must have both records, both suffixed
            self.assertIn('>Syllis_okadai__A', nt)
            self.assertIn('>Syllis_okadai__B', nt)
            # aa must also suffix sample A -- key regression assertion
            self.assertIn('>Syllis_okadai__A', aa,
                          'aa header must be suffixed to match nt; got: %r' % aa)
            self.assertNotIn('>Syllis_okadai\n', aa,
                             'bare header in aa would diverge from nt suffix')

    def test_no_collision_no_suffix(self):
        """Two samples with DISTINCT species headers must not gain any __ suffix."""
        with tempfile.TemporaryDirectory() as root:
            seqs = os.path.join(root, 'seqs')
            occ = os.path.join(root, 'occ')
            write(os.path.join(seqs, 'COX1__A.fna'), '>Syllis_okadai\nAAAA\n')
            write(os.path.join(seqs, 'COX1__A.faa'), '>Syllis_okadai\nMK\n')
            write(os.path.join(seqs, 'COX1__B.fna'), '>Syllis_gracilis\nCCCC\n')
            write(os.path.join(seqs, 'COX1__B.faa'), '>Syllis_gracilis\nML\n')
            write(os.path.join(occ, 'dummy.occupancy.tsv'),
                  'sample\tspecies\tgene\tlen_nt\tlen_aa\tflags\n')
            out = os.path.join(root, 'out')
            subprocess.run(
                [sys.executable, os.path.join(BIN, 'aggregate_mito_genes.py'),
                 '--seqs-dir', seqs, '--occ-dir', occ, '--outdir', out],
                check=True)
            with open(os.path.join(out, 'nt', 'COX1.fasta')) as fh:
                nt = fh.read()
            with open(os.path.join(out, 'aa', 'COX1.faa')) as fh:
                aa = fh.read()
            for line in (nt + aa).splitlines():
                if line.startswith('>'):
                    self.assertNotIn('__', line,
                                     'no suffix expected for distinct headers; got: %r' % line)


if __name__ == '__main__':
    unittest.main()
