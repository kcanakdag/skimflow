import os
import subprocess
import sys
import tempfile
import unittest

BIN = os.path.join(os.path.dirname(__file__), '..', 'bin')
FIX = os.path.join(os.path.dirname(__file__), 'fixtures', 'mitos2_OkadaiA')
sys.path.insert(0, BIN)
import harvest_mito_genes as h  # noqa: E402


class TestCanonical(unittest.TestCase):
    def test_synonyms_map_to_canonical(self):
        self.assertEqual(h.canonical_gene('cox1'), 'cox1')
        self.assertEqual(h.canonical_gene('COI'), 'cox1')
        self.assertEqual(h.canonical_gene('MT-CO1'), 'cox1')
        self.assertEqual(h.canonical_gene('cytb'), 'cob')
        self.assertEqual(h.canonical_gene('ND4L'), 'nad4l')
        self.assertEqual(h.canonical_gene('16S'), 'rrnl')
        self.assertEqual(h.canonical_gene('s-rRNA'), 'rrns')

    def test_nad4l_not_confused_with_nad4(self):
        self.assertEqual(h.canonical_gene('nad4l'), 'nad4l')
        self.assertEqual(h.canonical_gene('nad4'), 'nad4')

    def test_non_target_returns_none(self):
        self.assertIsNone(h.canonical_gene('trnM'))
        self.assertIsNone(h.canonical_gene('OH'))

    def test_split_tags(self):
        self.assertEqual(h.split_tags('cox1_1'), ('cox1', None))
        self.assertEqual(h.split_tags('cox1copy2'), ('cox1', None))
        self.assertEqual(h.split_tags('rrnL-a'), ('rrnL', 'a'))
        self.assertEqual(h.split_tags('nad4l'), ('nad4l', None))


class TestParseAndSelect(unittest.TestCase):
    def setUp(self):
        with open(os.path.join(FIX, 'result.fas')) as fh:
            self.nt = list(h.parse_mitos_fasta(fh.read()))
        with open(os.path.join(FIX, 'result.faa')) as fh:
            self.aa = list(h.parse_mitos_fasta(fh.read()))

    def test_parse_fields(self):
        first = self.nt[0]
        self.assertEqual(first['seqid'], 'Syllis_okadai')
        self.assertEqual((first['start'], first['end']), (1, 30))
        self.assertEqual(first['strand'], '+')
        self.assertEqual(first['canon'], 'cox1')

    def test_longest_copy_kept(self):
        feats = h.select_features(self.nt)
        self.assertEqual(len(feats['cox1']['seq']), 30)  # not the 16 bp cox1_1
        self.assertFalse(feats['cox1']['merged'])

    def test_rrna_fragments_merged(self):
        feats = h.select_features(self.nt)
        self.assertTrue(feats['rrnl']['merged'])
        # Minus strand: coding order is descending genomic coordinate, so the
        # higher-coord fragment (rrnL-a at 715-730) comes first, not letter order.
        self.assertEqual(feats['rrnl']['seq'], 'G' * 16 + 'A' * 13)

    def test_trna_excluded(self):
        feats = h.select_features(self.nt)
        self.assertNotIn(None, feats)
        self.assertEqual(set(feats), {'cox1', 'atp8', 'nad4l', 'rrnl'})

    def test_aa_coupled_to_same_locus(self):
        feats = h.select_features(self.nt)
        aa = h.select_aa(feats, self.aa)
        self.assertEqual(aa['cox1'], 'MAYGKF*')   # 30 bp locus, not the cox1_1 copy
        self.assertIn('nad4l', aa)

    def test_internal_stop_detection(self):
        self.assertTrue(h.has_internal_stop('MP*QLN'))
        self.assertFalse(h.has_internal_stop('MAYGKF*'))

    def test_sanitize_label(self):
        self.assertEqual(h.sanitize_label('Syllis okadai'), 'Syllis_okadai')
        # whitespace -> '_' happens before illegal-char stripping, so the space
        # after 'sp.' becomes '_' and the '(', ')', ';' are dropped.
        self.assertEqual(h.sanitize_label('sp. (A);1'), 'sp._A1')


class TestCli(unittest.TestCase):
    def test_end_to_end_writes_files(self):
        with tempfile.TemporaryDirectory() as out:
            subprocess.run(
                [sys.executable, os.path.join(BIN, 'harvest_mito_genes.py'),
                 '--mitos-dir', FIX, '--sample-id', 'OkadaiA',
                 '--species-id', 'Syllis okadai', '--outdir', out],
                check=True)
            seqs = os.path.join(out, 'seqs')
            self.assertTrue(os.path.exists(os.path.join(seqs, 'COX1__OkadaiA.fna')))
            self.assertTrue(os.path.exists(os.path.join(seqs, 'COX1__OkadaiA.faa')))
            self.assertTrue(os.path.exists(os.path.join(seqs, '16S__OkadaiA.fna')))
            self.assertFalse(os.path.exists(os.path.join(seqs, '16S__OkadaiA.faa')))  # rRNA has no aa
            with open(os.path.join(seqs, 'COX1__OkadaiA.fna')) as fh:
                self.assertEqual(fh.readline().strip(), '>Syllis_okadai')
            occ = os.path.join(out, 'OkadaiA.occupancy.tsv')
            self.assertTrue(os.path.exists(occ))
            with open(occ) as fh:
                text = fh.read()
            self.assertIn('merged_fragments', text)
            self.assertIn('internal_stop', text)


if __name__ == '__main__':
    unittest.main()
