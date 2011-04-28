"""Generates a csv file representing the ratio of leaves to edges in a tree."""

from Bio import Phylo

from collections import Counter, defaultdict
import argparse
import csv
import sys

from taxtastic import algotax, refpkg

def build_parser(parser):
    parser.add_argument('refpkg', nargs=1,
        help='the reference package to operate on')
    parser.add_argument('-o', '--outfile',
        type=argparse.FileType('w'), default='-',
        help='the name of the csv file to write out')

def action(args):
    args.outfile = csv.writer(args.outfile)
    args.outfile.writerow(['rank', 'name', 'leaves', 'edges'])
    rp = refpkg.Refpkg(args.refpkg[0])
    rp.load_db()
    with rp.resource('tree_file', 'rU') as fobj:
        tree = Phylo.read(fobj, 'newick')

    curs = rp.db.cursor()
    seq_colors = curs.execute("""
        SELECT t.rank,
               seqname,
               t.tax_name
        FROM   parents
               JOIN taxa t
                 ON t.tax_id = parent
               JOIN sequences s
                 ON s.tax_id = child
    """)
    rank_map = defaultdict(list)
    for rank, seq, color in seq_colors:
        rank_map[rank].append((seq, color))

    rank_order = curs.execute("""
        SELECT   rank
        FROM     ranks
        ORDER BY rank_order ASC
    """)

    for rank, in rank_order:
        seq_colors = rank_map[rank]
        if not seq_colors:
            continue
        clade_map = {c.name: c for c in tree.get_terminals()}
        colors = {}
        color_leaf_counts = Counter()
        for seq, color in seq_colors:
            colors[clade_map[seq]] = color
            color_leaf_counts[color] += 1
        metadata = algotax.color_clades(tree, colors)
        color_edge_counts = Counter()
        for edge_colors in metadata.cut_colors.itervalues():
            color_edge_counts += Counter(edge_colors)

        for color in (color_leaf_counts.viewkeys() |
                      color_edge_counts.viewkeys()):
            args.outfile.writerow([
                rank,
                color,
                color_leaf_counts[color],
                color_edge_counts[color]
            ])
