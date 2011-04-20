from Bio import Phylo
from Bio.Phylo.PhyloXML import BranchColor

from StringIO import StringIO
import collections
import functools
import itertools
import argparse
import operator
import os.path
import pprint
import copy
import sys

import algotax
import refpkg

def main(argv):
    parser = argparse.ArgumentParser(
        description='find discordance in a reference package')
    parser.add_argument('refpkg', nargs=1,
        help='the reference package to operate on')
    parser.add_argument('-o', '--outdir',
        help='if specified, write discordance trees in phyloxml format '
        'to this directory')
    parser.add_argument('-s', '--summary',
        action='store_const', dest='verbosity', const=1, default=2,
        help='show only a summary instead of the full output')

    args = parser.parse_args(argv)
    def out(level, msg='', *params):
        if args.verbosity >= level:
            print msg % params
    p = functools.partial(out, 2)

    rp = refpkg.Refpkg(args.refpkg[0])

    sio = StringIO()
    Phylo.convert(rp.resource('tree_file', 'rU'), 'newick', sio, 'phyloxml')
    sio.seek(0)
    tree = Phylo.read(sio, 'phyloxml')

    rp.load_db()
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
    rank_map = collections.defaultdict(list)
    for rank, seq, color in seq_colors:
        rank_map[rank].append((seq, color))

    rank_order = curs.execute("""
        SELECT   rank
        FROM     ranks
        ORDER BY rank_order ASC
    """)

    discordance_trees = []
    for rank, in rank_order:
        seq_colors = rank_map[rank]
        if not seq_colors:
            continue
        clade_map = {c.name: c for c in tree.get_terminals()}
        colors = {clade_map[seq]: color for seq, color in seq_colors}
        p('calculating for %s (|coloring| is %d)', rank, len(colors))
        metadata = algotax.color_clades(tree, colors)
        badness = sum(len(cut) - 1
            for node, cut in metadata.cut_colors.iteritems()
            if len(cut) > 1 and node != tree.root)
        if not badness:
            p('  completely convex; skipping')
            continue

        p('  badness: %d', badness)
        results = algotax.walk(tree.root, metadata)
        p('  possibilities:')
        counter = itertools.count()
        discordance_tree = None
        for k in sorted(results, key=results.get):
            p('    %d [%s]', len(results[k]), '; '.join(k) or '(none)')
            if not k:
                continue

            if args.outdir:
                # XXX find a better method of doing this
                discordance_tree = copy.deepcopy(tree)

            cut_nodes = set(tree.get_terminals()) - results[k]
            for node in sorted(cut_nodes, key=operator.attrgetter('name')):
                p('      %s', node.name)

                if not discordance_tree:
                    continue
                node = discordance_tree.find_any(name=node.name)
                node.color = BranchColor(255, 0, 0)
                node.width = 5

            if not discordance_tree:
                continue
            fname = 'discordance-%s-%d.xml' % (rank, next(counter))
            Phylo.write(
                discordance_tree, os.path.join(args.outdir, fname), 'phyloxml')
            discordance_trees.append(fname)

        p('    %d [whole set]', len(clade_map))
    if discordance_trees:
        if args.verbosity > 1:
            out(1)
        out(1, 'discordance trees written:')
        for d in discordance_trees:
            out(1, '  %s', d)

if __name__ == '__main__':
    main(sys.argv[1:])
