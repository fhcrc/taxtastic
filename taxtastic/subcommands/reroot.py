"""Reroots a reference package"""
import tempfile
import logging
import os

from Bio import Phylo

from taxtastic import algotax, refpkg

log = logging.getLogger(__name__)

def build_parser(parser):
    parser.add_argument('refpkg', nargs=1,
        help='the reference package to operate on')
    parser.add_argument('-p', '--pretend',
        action='store_true', default=False,
        help="don't save the rerooted tree; just attempt the rerooting.")

def action(args):
    log.info('loading reference package')
    rp = refpkg.Refpkg(args.refpkg[0])
    rp.load_db()
    with rp.resource('tree_file', 'rU') as fobj:
        tree = Phylo.read(fobj, 'newick')
    log.info('rerooting')
    root = algotax.reroot_from_rp(tree.root, rp)
    if log.isEnabledFor(logging.DEBUG):
        log.debug('new root is %d steps from the old root',
                  len(tree.get_path(root)))
    tree.root_with_outgroup(root)

    if args.pretend:
        return
    log.info('saving reference package')
    fd, name = tempfile.mkstemp()
    with os.fdopen(fd, 'w') as fobj:
        Phylo.write(tree, fobj, 'newick',
                    branchlengths_only=True,
                    format_branch_length='%0.6f')
    os.rename(name, rp.resource_path('tree_file'))
    rp.rehash('tree_file')
    rp.save()
