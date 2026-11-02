import os
import subprocess
import sys
import tempfile
import unittest

BIN = os.path.join(os.path.dirname(__file__), '..', 'bin')
SCRIPT = os.path.join(BIN, 'tree_to_svg.py')
sys.path.insert(0, BIN)
import tree_to_svg as t  # noqa: E402

# Real IQ-TREE output (pilot run, nt supermatrix): trifurcating basal node,
# branch lengths and 'SH-aLRT/UFBoot' support on the internal node.
PILOT = ('(Eurysyllis_tuberculata:0.5903534313,(Syllis_okadai:0.6353816251,'
         'Virchowia_sp:0.4644782576)91.3/87:0.1766510580,'
         'Trypanosyllis_sp:0.1887325416);')


class TestParse(unittest.TestCase):
    def test_trifurcating_root(self):
        root = t.parse_newick(PILOT)
        self.assertEqual(len(root.children), 3)
        self.assertEqual(t.count_tips(root), 4)
        inner = [c for c in root.children if c.children][0]
        self.assertEqual(inner.support, '91.3/87')
        self.assertAlmostEqual(inner.length, 0.1766510580)
        names = sorted(c.name for c in root.children if c.is_leaf())
        self.assertEqual(names, ['Eurysyllis_tuberculata', 'Trypanosyllis_sp'])

    def test_quoted_label_does_not_crash(self):
        # Quoted labels are not unquoted (see parse_newick), but a quoted name
        # without separators inside must still parse into the right shape.
        root = t.parse_newick("('Syllis_okadai':0.1,(B:0.2,C:0.3)95/100:0.05,D:0.4);")
        self.assertEqual(t.count_tips(root), 4)
        self.assertEqual(root.children[0].name, "'Syllis_okadai'")
        self.assertAlmostEqual(root.children[0].length, 0.1)

    def test_malformed_raises(self):
        with self.assertRaises(ValueError):
            t.parse_newick('(A:0.1,B:0.2')


class TestCli(unittest.TestCase):
    def _render(self, newick, *extra):
        with tempfile.TemporaryDirectory() as d:
            tree = os.path.join(d, 'nt_species.treefile')
            out = os.path.join(d, 'nt_tree_mqc.html')
            with open(tree, 'w') as fh:
                fh.write(newick + '\n')
            subprocess.run(
                [sys.executable, SCRIPT, '--treefile', tree, '--seqtype', 'nt',
                 '--out', out] + list(extra),
                check=True)
            with open(out) as fh:
                return fh.read()

    def test_renders_unrooted_caption_and_tips(self):
        text = self._render(PILOT)
        self.assertTrue(text.startswith('<!--'))
        self.assertIn('Unrooted maximum-likelihood species tree', text)
        self.assertIn('SH-aLRT / UFBoot support', text)
        self.assertIn('4 taxa.', text)
        self.assertIn('<svg', text)
        self.assertIn('Syllis okadai', text)   # '_' shown as a space
        self.assertIn('>91.3/87</text>', text)

    def test_support_label_none(self):
        text = self._render(PILOT, '--support-label', 'none')
        self.assertIn('no branch support was computed', text)


if __name__ == '__main__':
    unittest.main()
