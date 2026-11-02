import os
import sys
import tempfile
import unittest

BIN = os.path.join(os.path.dirname(__file__), '..', 'bin')
FIX = os.path.join(os.path.dirname(__file__), 'fixtures', 'metrics')
sys.path.insert(0, BIN)
import collect_metrics as cm  # noqa: E402

SUMMARY = (
    '# BUSCO version is: 5.x\n'
    '# Summarized benchmarking in BUSCO notation for file /w/{id}.contigs.fasta\n'
    '\tC:31.0%[S:30.0%,D:1.0%],F:2.0%,M:67.0%,n:954\n'
)


class TestBuscoId(unittest.TestCase):
    def test_dotted_id_kept_whole(self):
        with tempfile.TemporaryDirectory() as d:
            with open(os.path.join(d, 'short_summary.specific.metazoa_odb10.x_busco.txt'), 'w') as fh:
                fh.write(SUMMARY.format(id='Trypanosyllis_sp.1'))
            table = {}
            cm.parse_busco(d, table)
            self.assertEqual(list(table), ['Trypanosyllis_sp.1'])
            self.assertEqual(table['Trypanosyllis_sp.1']['busco_C'], '31.0')

    def test_filename_fallback_keeps_dots(self):
        with tempfile.TemporaryDirectory() as d:
            name = 'short_summary.specific.metazoa_odb10.Trypanosyllis_sp.1_busco.txt'
            with open(os.path.join(d, name), 'w') as fh:
                fh.write('\tC:31.0%[S:30.0%,D:1.0%],F:2.0%,M:67.0%,n:954\n')
            table = {}
            cm.parse_busco(d, table)
            self.assertEqual(list(table), ['Trypanosyllis_sp.1'])


class TestRealFixtures(unittest.TestCase):
    # Copied from real runs: a BUSCO 6 short_summary (pilot Syllis_okadai, work
    # path shortened) and RESPECT MultiQC tables with and without an estimate.
    def test_real_busco_summary(self):
        table = {}
        cm.parse_busco(os.path.join(FIX, 'busco'), table)
        self.assertEqual(list(table), ['Syllis_okadai'])
        r = table['Syllis_okadai']
        self.assertEqual((r['busco_C'], r['busco_S'], r['busco_D'],
                          r['busco_F'], r['busco_M']),
                         ('73.2', '72.9', '0.3', '21.7', '5.1'))

    def test_real_respect_values_and_na_passthrough(self):
        table = {}
        cm.parse_respect(os.path.join(FIX, 'respect'), table)
        self.assertEqual(sorted(table), ['Paraopisthosyllis_rufa', 'Pmisa'])
        self.assertEqual(table['Pmisa']['genome_size_bp'], '426769590')
        self.assertEqual(table['Pmisa']['coverage_x'], '0.35')
        self.assertEqual(table['Pmisa']['expected_size_bp'], 'NA')
        rufa = table['Paraopisthosyllis_rufa']
        self.assertEqual((rufa['genome_size_bp'], rufa['coverage_x']), ('NA', 'NA'))

    def test_meta_fills_na_expected_size_from_respect(self):
        table = {}
        cm.parse_respect(os.path.join(FIX, 'respect'), table)
        with tempfile.TemporaryDirectory() as d:
            with open(os.path.join(d, 'samples_meta.tsv'), 'w') as fh:
                fh.write('Pmisa\tPlatynereis massiliensis\t1200000000\n')
            cm.parse_meta(d, table)
        self.assertEqual(table['Pmisa']['expected_size_bp'], '1200000000')
        self.assertEqual(table['Pmisa']['species'], 'Platynereis_massiliensis')
        self.assertEqual(table['Pmisa']['genome_size_bp'], '426769590')


class TestMeta(unittest.TestCase):
    def _run(self, table, lines):
        with tempfile.TemporaryDirectory() as d:
            with open(os.path.join(d, 'samples_meta.tsv'), 'w') as fh:
                fh.write('\n'.join(lines) + '\n')
            cm.parse_meta(d, table)
        return table

    def test_fills_only_blank_fields(self):
        table = {
            'demo': {'sample': 'demo', 'expected_size_bp': 'NA'},
            'okadai': {'sample': 'okadai', 'species': 'From_occ', 'expected_size_bp': '5'},
        }
        self._run(table, ['demo\tBostrychus sinensis\t1000',
                          'okadai\tOther\t9',
                          'ghost\tGhost_sp\t1'])
        self.assertEqual(table['demo']['species'], 'Bostrychus_sinensis')
        self.assertEqual(table['demo']['expected_size_bp'], '1000')
        self.assertEqual(table['okadai']['species'], 'From_occ')
        self.assertEqual(table['okadai']['expected_size_bp'], '5')
        self.assertNotIn('ghost', table)

    def test_unknown_species_skipped(self):
        table = {'s': {'sample': 's'}}
        self._run(table, ['s\tunknown\t'])
        self.assertNotIn('species', table['s'])
        self.assertNotIn('expected_size_bp', table['s'])


if __name__ == '__main__':
    unittest.main()
