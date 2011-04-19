from Bio import Phylo
from Bio.Phylo.PhyloXML import BranchColor

from StringIO import StringIO
import collections
import itertools
import operator
import algotax
import pprint
import refpkg
import copy
import sys

def main(refpkg_path):
    rp = refpkg.Refpkg(refpkg_path)

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

    for rank, in rank_order:
        seq_colors = rank_map[rank]
        if not seq_colors:
            continue
        clade_map = {c.name: c for c in tree.get_terminals()}
        colors = {clade_map[seq]: color for seq, color in seq_colors}
        print 'calculating for %s (|coloring| is %d)' % (rank, len(colors))
        metadata = algotax.color_clades(tree, colors)
        badness = sum(len(cut) - 1
            for node, cut in metadata.cut_colors.iteritems()
            if len(cut) > 1 and node != tree.root)
        if not badness:
            print '  completely convex; skipping'
            continue

        print '  badness: %d' % (badness,)
        results = algotax.walk(tree.root, metadata)
        print '  possibilities:'
        counter = itertools.count()
        for k in sorted(results, key=results.get):
            print '    %d [%s]' % (len(results[k]), '; '.join(k) or '(none)')
            if not k:
                continue
            # XXX find a better method of doing this
            discordance_tree = copy.deepcopy(tree)
            uncut_nodes = set(n.name for n in results[k])
            for node in sorted(discordance_tree.get_terminals(),
                    key=operator.attrgetter('name')):
                if node.name in uncut_nodes:
                    continue
                print '      %s' % (node.name)
                node.color = BranchColor(255, 0, 0)
                node.width = 5
            fname = 'discordance-%s-%d.xml' % (rank, next(counter))
            Phylo.write(discordance_tree, fname, 'phyloxml')
            print '      [discordance tree written to %s]' % (fname,)
        print '    %d [whole set]' % (len(clade_map),)

if __name__ == '__main__':
    main(*sys.argv[1:])
