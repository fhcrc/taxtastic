import matplotlib.pyplot as plt
from Bio import Phylo
import numpy as np

from collections import Counter, defaultdict
import argparse
import sys

import algotax
import refpkg

def main(argv):
    parser = argparse.ArgumentParser(
        description='reroot a reference package')
    parser.add_argument('refpkg', nargs=1,
        help='the reference package to operate on')

    args = parser.parse_args(argv)

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
        counts = [(color, color_leaf_counts[color], color_edge_counts[color])
            for color in
            color_leaf_counts.viewkeys() | color_edge_counts.viewkeys()]
        plot_color, plot_leaves, plot_edges = zip(*counts)
        plt.scatter(plot_leaves, plot_edges, marker='o', alpha=.6)
        ylim = max(max(plot_leaves) * 2, max(plot_edges)) * 1.1
        xlim = ylim / 2
        t = np.arange(0, xlim, .01)
        slope = 2 * t - 2
        plt.plot(t, slope)
        plt.xlabel('leaves')
        plt.ylabel('edges')
        plt.axis([1, xlim, 0, ylim])
        plt.savefig('%s.png' % rank.encode())
        plt.clf()

if __name__ == '__main__':
    main(sys.argv[1:])
