from StringIO import StringIO
from Bio import Phylo
import collections
import algotax
import pprint

def intersection(it):
    ret = None
    for x in it:
        if ret is None:
            ret = set(x)
        else:
            ret.intersection_update(x)
    return ret or set()

CladeMetadata = collections.namedtuple(
    'CladeMetadata', 'parents colors cut_colors')

def color_clades(tree, colors):
    parents = {tree.root: None}
    cut_colors = collections.defaultdict(set)
    stack = [('down', tree.root, None)]
    while stack:
        phase, cur, color = stack.pop()
        if phase == 'down':
            if not cur.clades:
                stack.append(('up', cur, cur.name))
            else:
                for child in cur.clades:
                    parents[child] = cur
                    stack.append(('down', child, None))
        elif phase == 'up':
            if cur is None or color in cut_colors[cur]:
                continue
            cut_colors[cur].add(color)
            stack.append(('up', parents[cur], color))

    # This isn't actually generalized to >2 children.
    # It's just shorter this way.
    stack = [(tree.root, set())]
    while stack:
        node, okayed = stack.pop()
        if not node.clades:
            continue
        okayed = intersection(cut_colors[e] for e in node.clades) | okayed
        for e in node.clades:
            e_ = cut_colors[e] & okayed
            if e_ != cut_colors[e]:
                stack.append((e, okayed))
                cut_colors[e] = e_
    return CladeMetadata(parents, colors, cut_colors)

def main():
    tree_string = StringIO('((A,B),(A,(B,B)))')
    tree = Phylo.read(tree_string, 'newick')
    colors = {n: n.name for n in tree.get_terminals()}
    print tree
    metadata = color_clades(tree, colors)
    print metadata
    print algotax.walk(tree.root, metadata)

if __name__ == '__main__':
    main()
