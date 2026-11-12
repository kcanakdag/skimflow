import os
import sys
import unittest

BENCH = os.path.join(os.path.dirname(__file__), '..', 'benchmarks')
sys.path.insert(0, BENCH)
import validate_mitogenomes as v  # noqa: E402


class TestNormGene(unittest.TestCase):
    def test_names(self):
        cases = {
            'ND4L': 'nad4l', 'NAD1': 'nad1', 'CYTB': 'cob', 'COX1': 'cox1',
            'cox1-a': 'cox1', 'trnL1(tag)': 'trnL', 'trnS2': 'trnS',
            'tRNA-Leu': 'trnL', '16S ribosomal RNA': 'rrnL', 'l-rRNA': 'rrnL',
            's-rRNA': 'rrnS', 'cytochrome c oxidase subunit I': 'cox1',
            'NADH dehydrogenase subunit 4L': 'nad4l',
            'ATP synthase F0 subunit 8': 'atp8', 'OH_1-a': None, 'gpII': None,
        }
        for raw, want in cases.items():
            self.assertEqual(v.norm_gene(raw), want, raw)

    def test_expected_set_is_37(self):
        self.assertEqual(sum(v.expected_counts().values()), 37)


class TestFeatureTable(unittest.TestCase):
    FT = ('>Feature gb|X|\n'
          '10\t80\tgene\n\t\t\tgene\ttrnL2\n'
          '10\t80\ttRNA\n\t\t\tproduct\ttRNA-Leu\n'
          '100\t400\tCDS\n'          # joined CDS without a gene feature
          '900\t1200\n'
          '\t\t\tproduct\tcytochrome c oxidase subunit I\n'
          '2000\t1500\tgene\n\t\t\tgene\tND6\n'
          '2000\t1500\tCDS\n\t\t\tproduct\tNADH dehydrogenase subunit 6\n')

    def test_parse(self):
        feats = v.parse_feature_table(self.FT)
        self.assertEqual(feats, [(10, 80, '+', 'trnL'), (100, 1200, '+', 'cox1'),
                                 (1500, 2000, '-', 'nad6')])

    # skim2mito/ENA style: genes carry only a locus_tag, protein genes are
    # mRNA features, the name sits in a note "ID:transcript_<gene>"
    FT_ENA = ('>Feature emb|OZ1|\n'
              '1\t70\tgene\n\t\t\tlocus_tag\tL1\n'
              '1\t70\ttRNA\n\t\t\tnote\tID:transcript_trnL1;~source:mitfi\n'
              '1\t70\texon\n\t\t\tstandard_name\ttrnL1\n'
              '100\t900\tgene\n\t\t\tlocus_tag\tL2\n'
              '100\t900\tmRNA\n\t\t\tnote\tID:transcript_nad1;~source:mitos\n'
              '1000\t1900\trRNA\n\t\t\tnote\tID:transcript_rrnL;~source:mitfi\n')

    def test_parse_ena_style(self):
        feats = v.parse_feature_table(self.FT_ENA)
        self.assertEqual(feats, [(1, 70, '+', 'trnL'), (100, 900, '+', 'nad1'),
                                 (1000, 1900, '+', 'rrnL')])

    def test_mrna_ignored_when_cds_present(self):
        ft = ('>Feature gb|X|\n'
              '100\t900\tmRNA\n\t\t\tproduct\tNADH dehydrogenase subunit 1\n'
              '100\t900\tCDS\n\t\t\tproduct\tNADH dehydrogenase subunit 1\n')
        self.assertEqual(v.parse_feature_table(ft), [(100, 900, '+', 'nad1')])


class TestOrder(unittest.TestCase):
    def test_reverse_reading_is_same_adjacency(self):
        fwd = [('cox1', '+'), ('cox2', '+'), ('nad6', '-')]
        rev = [('nad6', '+'), ('cox2', '-'), ('cox1', '-')]
        self.assertEqual(v.adjacencies(fwd, True), v.adjacencies(rev, True))

    def test_linear_loses_wrap(self):
        o = [('a', '+'), ('b', '+'), ('c', '+')]
        self.assertEqual(len(v.adjacencies(o, True)), 3)
        self.assertEqual(len(v.adjacencies(o, False)), 2)

    def test_split_parts_collapse(self):
        feats = [(1, 5, '+', 'cox1'), (9, 12, '+', 'cox1'), (20, 30, '+', 'cox2')]
        self.assertEqual(v.order_of(feats), [('cox1', '+'), ('cox2', '+')])

    def test_adjacent_trna_pair_not_collapsed(self):
        # gastropods: nad1 trnL2 trnL1 rrnL; both leucine tRNAs must count
        feats = [(1, 900, '-', 'nad1'), (901, 970, '-', 'trnL'),
                 (971, 1040, '-', 'trnL'), (1041, 2400, '-', 'rrnL')]
        self.assertEqual(v.gene_counts(v.collapse_parts(feats))['trnL'], 2)


class TestBlastStats(unittest.TestCase):
    def test_coverage_per_sequence(self):
        # two contigs, each fully aligned, same query coordinates
        hsps = [(100.0, 100, 1, 100, 1, 100, 200.0, 'c1'),
                (99.0, 100, 1, 100, 201, 300, 190.0, 'c2')]
        ref_cov, asm_cov, ident = v.blast_stats(hsps, qlen=200, slen=400)
        self.assertAlmostEqual(asm_cov, 1.0)
        self.assertAlmostEqual(ref_cov, 0.5)
        self.assertAlmostEqual(ident, 99.5)

    def test_repeats_not_counted_twice_for_identity(self):
        hsps = [(100.0, 300, 1, 300, 1, 300, 600.0, 'c1'),
                (90.0, 300, 301, 600, 1, 300, 400.0, 'c1')]
        _, _, ident = v.blast_stats(hsps, qlen=600, slen=300)
        self.assertAlmostEqual(ident, 100.0)


if __name__ == '__main__':
    unittest.main()
