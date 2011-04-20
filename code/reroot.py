from Bio import Phylo

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
    root = algotax.reroot(tree.root, rp)
    tree.root_with_outgroup(root)
    with rp.resource('tree_file', 'w') as fobj:
        Phylo.write(tree, fobj, 'newick')
    rp.rehash('tree_file')
    rp.save()

if __name__ == '__main__':
    main(sys.argv[1:])
