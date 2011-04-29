"""Generates a csv file representing the ratio of leaves to edges in a tree."""

from Bio import Phylo

from collections import Counter, defaultdict
from itertools import combinations
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
    args.outfile.writerow(['rank', 'name', 'leaves', 'mrca_leaves'])
    rp = refpkg.Refpkg(args.refpkg[0])
    rp.load_db()
    with rp.resource('tree_file', 'rU') as fobj:
        tree = Phylo.read(fobj, 'newick')
    clade_map = {c.name: c for c in tree.get_terminals()}

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
        colors = {}
        color_leaf_counts = Counter()
        for seq, color in seq_colors:
            colors[clade_map[seq]] = color
            color_leaf_counts[color] += 1
        metadata = algotax.color_clades(tree, colors)
        mrcas = {}
        for clade in tree.find_elements(order='level'):
            if clade is tree:
                continue
            colors_crossing = algotax.union(
                metadata.cut_colors[a] & metadata.cut_colors[b]
                for a, b in combinations(clade.clades, 2))
            for color in colors_crossing:
                mrcas.setdefault(color, clade)

        for color in color_leaf_counts.viewkeys() & mrcas.viewkeys():
            args.outfile.writerow([
                rank,
                color,
                color_leaf_counts[color],
                mrcas[color].count_terminals(),
            ])
