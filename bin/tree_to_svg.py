#!/usr/bin/env python3
"""Render an IQ-TREE Newick .treefile as a MultiQC custom-content section.

Emits `<seqtype>_tree_mqc.html`: a leading MultiQC comment header plus an inline
SVG phylogram (branch lengths to scale, tip labels, internal-node support). Pure
standard library only, so it runs in the plain python:3.12 gene_container with
nothing to install and no network. MultiQC picks up any `*_mqc.html` file and
renders it as a section, reading id/section_name/description from the header.
"""
import argparse
import html
import math
import sys


class Node:
    __slots__ = ("name", "length", "support", "children", "x", "y", "tipindex")

    def __init__(self):
        self.name = None
        self.length = 0.0
        self.support = None
        self.children = []
        self.x = 0.0
        self.y = 0.0
        self.tipindex = 0

    def is_leaf(self):
        return not self.children


def parse_newick(text):
    """Minimal recursive-descent Newick parser.

    Handles branch lengths (:x), internal-node labels (IQ-TREE writes support as
    'SH-aLRT/UFBoot', e.g. 91.3/87), and multifurcating roots (unrooted trees).
    Quoted labels are not expected from IQ-TREE species trees and are not
    handled.
    """
    text = text.strip()
    n = len(text)
    i = 0

    def read_token():
        nonlocal i
        start = i
        while i < n and text[i] not in ":,();":
            i += 1
        return text[start:i]

    def read_number():
        nonlocal i
        start = i
        while i < n and text[i] not in ",();":
            i += 1
        try:
            return float(text[start:i])
        except ValueError:
            return 0.0

    def parse_clade():
        nonlocal i
        node = Node()
        if i < n and text[i] == "(":
            i += 1  # consume '('
            while True:
                node.children.append(parse_clade())
                if i < n and text[i] == ",":
                    i += 1
                    continue
                if i < n and text[i] == ")":
                    i += 1
                    break
                raise ValueError("malformed newick near position %d" % i)
            label = read_token()  # internal label = support
            if label:
                node.support = label
        else:
            node.name = read_token()
        if i < n and text[i] == ":":
            i += 1
            node.length = read_number()
        return node

    root = parse_clade()
    return root


def count_tips(node):
    if node.is_leaf():
        return 1
    return sum(count_tips(c) for c in node.children)


def ladderize(node):
    """Order children by clade size for the classic staircase look."""
    if node.is_leaf():
        return
    for c in node.children:
        ladderize(c)
    node.children.sort(key=count_tips)


def layout(root):
    tips = []

    def collect(node):
        if node.is_leaf():
            node.tipindex = len(tips)
            tips.append(node)
        else:
            for c in node.children:
                collect(c)

    collect(root)

    def set_y(node):
        if node.is_leaf():
            node.y = float(node.tipindex)
        else:
            for c in node.children:
                set_y(c)
            node.y = sum(c.y for c in node.children) / len(node.children)

    set_y(root)

    root.length = 0.0

    def set_x(node, acc):
        node.x = acc + (node.length or 0.0)
        for c in node.children:
            set_x(c, node.x)

    set_x(root, 0.0)
    return tips


def nice_scale(maxx):
    if maxx <= 0:
        return 0.1
    target = maxx / 4.0
    exp = math.floor(math.log10(target))
    base = 10 ** exp
    for m in (1, 2, 5, 10):
        if m * base >= target:
            return m * base
    return 10 * base


def max_x(node):
    m = node.x
    for c in node.children:
        m = max(m, max_x(c))
    return m


def render_svg(root, tips):
    LINE = "#34495e"
    TXT = "#2c3e50"
    SUP = "#7f8c8d"
    row_h = 34
    top = 26
    bottom = 40
    left = 12
    branch_w = 430          # horizontal pixels for the branches themselves
    label_w = 250           # reserved space for tip labels on the right
    mx = max_x(root)
    x_scale = branch_w / mx if mx > 0 else 1.0
    height = top + max(1, len(tips)) * row_h + bottom
    width = left + branch_w + label_w

    def px(node):
        return left + node.x * x_scale

    def py(node):
        return top + node.y * row_h

    parts = []
    parts.append(
        '<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 %d %d" '
        'style="max-width:%dpx;width:100%%;height:auto;font-family:'
        '-apple-system,Segoe UI,Helvetica,Arial,sans-serif" '
        'role="img" aria-label="phylogenetic tree">' % (width, height, width)
    )

    # Edges: for every internal node draw a vertical connector at its x spanning
    # its children, and a horizontal branch from the node to each child.
    def draw_edges(node):
        if node.children:
            ys = [py(c) for c in node.children]
            parts.append(
                '<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" stroke="%s" '
                'stroke-width="1.6"/>' % (px(node), min(ys), px(node), max(ys), LINE)
            )
            for c in node.children:
                parts.append(
                    '<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" stroke="%s" '
                    'stroke-width="1.6"/>' % (px(node), py(c), px(c), py(c), LINE)
                )
                draw_edges(c)

    draw_edges(root)

    # Internal-node support labels (skip the root).
    def draw_support(node, is_root):
        if node.children:
            if not is_root and node.support:
                parts.append(
                    '<text x="%.1f" y="%.1f" font-size="10" fill="%s" '
                    'text-anchor="end">%s</text>'
                    % (px(node) - 3, py(node) - 4, SUP, html.escape(node.support))
                )
            for c in node.children:
                draw_support(c, False)

    draw_support(root, True)

    # Tip labels (italic species names).
    for t in tips:
        parts.append(
            '<text x="%.1f" y="%.1f" font-size="13" font-style="italic" '
            'fill="%s" dominant-baseline="middle">%s</text>'
            % (px(t) + 6, py(t), TXT, html.escape(t.name or ""))
        )

    # Scale bar (substitutions/site).
    sv = nice_scale(mx)
    bar = sv * x_scale
    by = height - bottom + 22
    bx = left
    parts.append(
        '<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" stroke="%s" '
        'stroke-width="1.6"/>' % (bx, by, bx + bar, by, TXT)
    )
    parts.append(
        '<text x="%.1f" y="%.1f" font-size="11" fill="%s" text-anchor="middle">'
        "%s</text>" % (bx + bar / 2, by + 14, TXT, _fmt(sv))
    )
    parts.append(
        '<text x="%.1f" y="%.1f" font-size="10" fill="%s">substitutions/site</text>'
        % (bx + bar + 8, by + 4, SUP)
    )
    parts.append("</svg>")
    return "".join(parts)


def _fmt(v):
    s = ("%.4f" % v).rstrip("0").rstrip(".")
    return s


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--treefile", required=True)
    ap.add_argument("--seqtype", required=True, help="nt or aa")
    ap.add_argument("--out", required=True, help="output <seqtype>_tree_mqc.html")
    args = ap.parse_args()

    with open(args.treefile) as fh:
        newick = fh.read().strip()

    if not newick:
        sys.exit("empty treefile: %s" % args.treefile)

    root = parse_newick(newick)
    ladderize(root)
    tips = layout(root)

    st = args.seqtype
    seqlabel = {"nt": "nucleotide", "aa": "amino-acid"}.get(st, st)
    svg = render_svg(root, tips)

    header = (
        "<!--\n"
        "id: 'phylogeny_%s_tree'\n"
        "section_name: 'Phylogeny tree (%s)'\n"
        "description: 'Maximum-likelihood species tree from the concatenated %s "
        "(%s) supermatrix (IQ-TREE, model-selected). Branch lengths are "
        "substitutions/site; internal labels are SH-aLRT / UFBoot support. "
        "%d taxa.'\n"
        "-->\n" % (st, st, st, seqlabel, len(tips))
    )
    with open(args.out, "w") as fh:
        fh.write(header)
        fh.write('<div style="padding:6px 2px 2px">')
        fh.write(svg)
        fh.write("</div>\n")


if __name__ == "__main__":
    main()
