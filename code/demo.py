from StringIO import StringIO
from Bio import Phylo
import algotax

def main():
    tree_string = StringIO('((A,B),(A,(B,B)))')
    tree = Phylo.read(tree_string, 'newick')
    colors = {n: n.name for n in tree.get_terminals()}
    print tree
    metadata = algotax.color_clades(tree, colors)
    print metadata
    print algotax.walk(tree.root, metadata)

if __name__ == '__main__':
    main()
