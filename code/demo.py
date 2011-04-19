from Bio import Phylo
import collections
import algotax
import pprint
import refpkg
import sys

def main(refpkg_path):
    rp = refpkg.Refpkg(refpkg_path)
    tree = Phylo.read(rp.resource('tree_file', 'rU'), 'newick')
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
        for k in sorted(results, key=results.get):
            print '    %d [%s]' % (len(results[k]), '; '.join(k) or '(none)')
            if not k:
                continue
            cut_nodes = set(tree.get_terminals()) - results[k]
            for nodename in sorted(node.name for node in cut_nodes):
                print '      %s' % (nodename)
        print '    %d [whole set]' % (len(clade_map),)

if __name__ == '__main__':
    main(*sys.argv[1:])
