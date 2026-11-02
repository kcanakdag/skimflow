import os
import subprocess
import sys
import tempfile
import unittest

BIN = os.path.join(os.path.dirname(__file__), '..', 'bin')
FIX = os.path.join(os.path.dirname(__file__), 'fixtures', 'mitos2_OkadaiA')
sys.path.insert(0, BIN)
import harvest_mito_genes as h  # noqa: E402

# Ten records copied verbatim from a real MITOS2 result.fas (D. melanogaster
# test run): OH_N and OH_N-x control-region tags, trnX(anticodon) names, and
# two PCGs.
REAL_FIX = os.path.join(os.path.dirname(__file__), 'fixtures', 'mitos2_dmel_real')


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
        self.assertEqual(h.split_tags('cox1_1'), ('cox1', None, '1'))
        self.assertEqual(h.split_tags('cox1copy2'), ('cox1', None, '2'))
        self.assertEqual(h.split_tags('rrnL-a'), ('rrnL', 'a', None))
        self.assertEqual(h.split_tags('nad4l'), ('nad4l', None, None))
        self.assertEqual(h.split_tags('rrnL_1-a'), ('rrnL', 'a', '1'))

    def test_copy_and_fragment_tags_combined(self):
        # MITOS2 names a feature NAME_COPY-PART.
        self.assertEqual(h.canonical_gene('rrnL_1-a'), 'rrnl')
        self.assertEqual(h.canonical_gene('rrnL_0-b'), 'rrnl')
        self.assertEqual(h.canonical_gene('cox1_1-a'), 'cox1')
        self.assertEqual(h.canonical_gene('nad4l_1'), 'nad4l')
        self.assertIsNone(h.canonical_gene('OH_2-a'))


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


class TestRealMitos2Headers(unittest.TestCase):
    def setUp(self):
        with open(os.path.join(REAL_FIX, 'result.fas')) as fh:
            self.nt = list(h.parse_mitos_fasta(fh.read()))

    def test_all_records_parsed(self):
        self.assertEqual(len(self.nt), 10)
        for r in self.nt:
            self.assertEqual(r['seqid'], 'scaffold_1--1676+')
            self.assertEqual(len(r['seq']), r['end'] - r['start'] + 1)

    def test_oh_and_trna_map_to_none(self):
        for r in self.nt:
            if r['raw'].startswith(('OH', 'trn')):
                self.assertIsNone(r['canon'], r['raw'])

    def test_oh_combined_tag_split(self):
        oh = {r['raw']: r for r in self.nt if r['raw'].startswith('OH')}
        self.assertEqual((oh['OH_2-c']['frag'], oh['OH_2-c']['copy']), ('c', '2'))
        self.assertEqual((oh['OH_1']['frag'], oh['OH_1']['copy']), (None, '1'))

    def test_only_target_genes_selected(self):
        feats = h.select_features(self.nt)
        self.assertEqual(set(feats), {'nad4l', 'atp8'})
        self.assertEqual((feats['nad4l']['start'], feats['nad4l']['end']), (6840, 7130))
        self.assertEqual(len(feats['atp8']['seq']), 162)
        self.assertFalse(feats['atp8']['merged'])


def _rec(raw, start, end, strand, seq):
    base, frag, copy = h.split_tags(raw)
    return {'seqid': 'c1', 'start': start, 'end': end, 'strand': strand,
            'raw': raw, 'canon': h.canonical_gene(raw), 'frag': frag,
            'copy': copy, 'seq': seq}


class TestCopyFragmentSelection(unittest.TestCase):
    def test_full_copy_beats_shorter_fragmented_copy(self):
        recs = [_rec('rrnL_0', 1, 100, '+', 'A' * 100),
                _rec('rrnL_1-a', 200, 229, '+', 'C' * 30),
                _rec('rrnL_1-b', 240, 269, '+', 'G' * 30)]
        feats = h.select_features(recs)
        self.assertEqual(feats['rrnl']['seq'], 'A' * 100)
        self.assertFalse(feats['rrnl']['merged'])

    def test_fragments_of_different_copies_not_joined(self):
        recs = [_rec('rrnS_0-a', 10, 49, '+', 'A' * 40),
                _rec('rrnS_0-b', 60, 99, '+', 'C' * 40),
                _rec('rrnS_1-a', 300, 349, '+', 'G' * 50),
                _rec('rrnS_1-b', 400, 439, '+', 'T' * 40)]
        feats = h.select_features(recs)
        self.assertEqual(feats['rrns']['seq'], 'G' * 50 + 'T' * 40)
        self.assertTrue(feats['rrns']['merged'])
        self.assertEqual((feats['rrns']['start'], feats['rrns']['end']), (300, 439))

    def test_combined_tag_gene_not_reported_absent(self):
        recs = [_rec('nad5_1-a', 10, 39, '+', 'A' * 30),
                _rec('nad5_1-b', 50, 79, '+', 'C' * 30)]
        feats = h.select_features(recs)
        self.assertEqual(feats['nad5']['seq'], 'A' * 30 + 'C' * 30)


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
            # Verify occupancy TSV structure
            lines = text.splitlines()
            self.assertEqual(lines[0], 'sample\tspecies\tgene\tlen_nt\tlen_aa\tflags')
            self.assertEqual(len(lines), 16)  # 1 header + 15 genes
            expected_genes = {'COX1', 'COX2', 'COX3', 'COB', 'NAD1', 'NAD2', 'NAD3', 'NAD4',
                              'NAD4L', 'NAD5', 'NAD6', 'ATP6', 'ATP8', '16S', '12S'}
            for line in lines[1:]:
                fields = line.split('\t')
                self.assertIn(fields[2], expected_genes)


if __name__ == '__main__':
    unittest.main()
