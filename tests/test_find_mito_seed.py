import os
import sys
import tempfile
import unittest

BIN = os.path.join(os.path.dirname(__file__), '..', 'bin')
sys.path.insert(0, BIN)
import find_mito_seed as f  # noqa: E402


def hit(gene, contig, qstart=1, qend=500, qlen=510, bits=800.0, sseq='MLWF'):
    return '\t'.join(map(str, [f'{gene}|Sp_x|P1', contig, 90, qend - qstart + 1,
                               qstart, qend, 1, 1500, 1e-50, bits, qlen, sseq]))


CONTIGS = ('>k141_1 flag=1 multi=4000.0 len=6000\n' + 'A' * 6000 + '\n'
           '>k141_2 flag=1 multi=1500.0 len=3000\n' + 'C' * 3000 + '\n'
           '>k141_3 flag=1 multi=90.0 len=1600\n' + 'G' * 1600 + '\n'
           '>k141_4 flag=1 multi=2.0 len=1500\n' + 'T' * 1500 + '\n'
           '>k141_5 flag=1 multi=1.0 len=400\n' + 'A' * 400 + '\n')


class TestFindMitoSeed(unittest.TestCase):
    def run_case(self, hit_lines, contigs=CONTIGS, extra=()):
        with tempfile.TemporaryDirectory() as d:
            h, c = os.path.join(d, 'hits.tsv'), os.path.join(d, 'ctg.fa')
            fa, tsv = os.path.join(d, 'seed.fa'), os.path.join(d, 's.tsv')
            open(h, 'w').write('\n'.join(hit_lines) + '\n')
            open(c, 'w').write(contigs)
            f.main(['--hits', h, '--contigs', c, '--sample', 'S1',
                    '--out-fasta', fa, '--out-tsv', tsv, *extra])
            seed = open(fa).read() if os.path.exists(fa) else None
            row = open(tsv).read().splitlines()[-1].split('\t')
        return seed, row

    def test_keeps_deep_clean_contigs_drops_contaminant_and_numt(self):
        seed, row = self.run_case([
            hit('cox1', 'k141_3', bits=900.0),           # contaminant: best score, low depth
            hit('cox1', 'k141_1', bits=850.0),           # own mito
            hit('nad5', 'k141_2', bits=600.0),           # own mito, other gene
            hit('cox1', 'k141_4', bits=700.0, sseq='ML*WF'),  # NUMT: stop codon
            hit('cob', 'k141_5', qstart=1, qend=50),     # too partial
        ])
        names = [l[1:] for l in seed.splitlines() if l.startswith('>')]
        self.assertEqual(names, ['k141_1', 'k141_2'])
        self.assertEqual(row[1:6], ['found', '2', '9000', '2/13', '4000.0'])
        self.assertEqual(row[8:10], ['1', '1'])
        self.assertIn('k141_3', row[10])

    def test_drops_contig_with_only_duplicate_genes(self):
        seed, row = self.run_case([
            hit('cox1', 'k141_1'), hit('nad1', 'k141_1'),
            hit('nad1', 'k141_2', bits=300.0),           # second copy of nad1
        ])
        names = [l[1:] for l in seed.splitlines() if l.startswith('>')]
        self.assertEqual(names, ['k141_1'])
        self.assertIn('duplicate genes', row[10])
        self.assertIn('k141_2', row[10])

    def test_none_found(self):
        seed, row = self.run_case([hit('cox1', 'k141_4', sseq='M*')])
        self.assertIsNone(seed)
        self.assertEqual(row[1], 'none')

    def test_max_seed_bp_caps_but_keeps_first(self):
        seed, row = self.run_case([hit('cox1', 'k141_1'), hit('cox2', 'k141_2')],
                                  extra=['--max-seed-bp', '5000'])
        self.assertEqual(row[2:4], ['1', '6000'])

    def test_without_depth_ranks_by_score(self):
        contigs = '>a\n' + 'A' * 100 + '\n>b\n' + 'C' * 100 + '\n'
        seed, row = self.run_case([hit('cox2', 'a', bits=10.0), hit('cox1', 'b', bits=99.0)],
                                  contigs=contigs)
        names = [l[1:] for l in seed.splitlines() if l.startswith('>')]
        self.assertEqual(names, ['b', 'a'])
        self.assertIn('no depth', row[10])


if __name__ == '__main__':
    unittest.main()
