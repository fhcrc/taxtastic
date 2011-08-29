"""Creates a reference package"""

import logging
import os
import time
import shutil
import hashlib
import re
import json
from collections import defaultdict

from taxtastic import utils
from taxtastic import package

log = logging.getLogger(__name__)

class ConfigError(Exception):
    pass

def build_parser(parser):
    parser.add_argument("-a", "--author",
        action="store", dest="author",
        help='Person who created the reference package', metavar='NAME')

    parser.add_argument("-c", "--clobber",
        action="store_true", dest="clobber", default = False,
        help= 'Delete an existing reference package.')

    parser.add_argument("-d", "--description",
        action="store", dest="description",
        help='An arbitrary description field', metavar='TEXT')

    parser.add_argument("-f", "--aln-fasta",
        action="store", dest="aln_fasta",
        help='Multiple alignment in fasta format', metavar='FILE')

    parser.add_argument("-i", "--seq-info",
        action="store", dest="seq_info",
        help='CSV format file describing the aligned reference ' + \
             'sequences, minimally containing the fields "seqname" ' + \
             'and "tax_id"', metavar='FILE')

    parser.add_argument("-l", "--locus",
        action="store", dest="locus", required=True,
        help='The locus described by the reference package', metavar='LOCUS')

    parser.add_argument("-m", "--mask",
        action="store", dest="mask",
        help='Text file containing a mask', metavar='FILE')

    parser.add_argument("-p", "--profile",
        action="store", dest="profile",
        help='Alignment profile', metavar='FILE')

    parser.add_argument('-P', '--package-name',
        action='store', dest='package_name',
        default='./taxtastic.refpkg', metavar='PATH',
        help='Name of output directory [default %(default)s]')

    parser.add_argument("-R", "--readme",
                        action="store", dest="readme",
                        help="README file describing the reference package",
                        metavar="FILE")

    parser.add_argument("-r", "--package-version",
        action="store", dest="package_version",
        help='Release version for the reference package', metavar='VERSION')

    parser.add_argument("-s", "--tree-stats",
        action="store", dest="tree_stats",
        help='File containing tree statistics (for example ' + \
             'RAxML_info.whatever")', metavar='FILE')

    parser.add_argument("-S", "--aln-sto",
        action="store", dest="aln_sto",
        help='Multiple alignment in Stockholm format', metavar='FILE')

    parser.add_argument("-t", "--tree-file",
        action="store", dest="tree",
        help='Phylogenetic tree in newick format',
        metavar='FILE')

    parser.add_argument("-T", "--taxonomy",
        action="store", dest="taxonomy",
        help='CSV format file defining the taxonomy. Fields include ' + \
             '"tax_id","parent_id","rank","tax_name" followed by a column ' + \
             'defining tax_id at each rank starting with root', metavar='FILE')


def action(args):
    """
    Create the reference package (a directory with a manifest named
    `manifest_name`).

     * args - output of argparse.ArgumentParser.parse_args (or
       presumably any other object with the required attributes:
       [TODO: list required attributes here])
     * MANIFEST_NAME - name of the JSON-format manifest file. Uses
       taxonomy.package.MANIFEST_NAME by default.
     * PACKAGE_CONTENTS - A dict defining sections and contents of
       each. Uses taxonomy.package.PACKAGE_CONTENTS by default.
     * PHYLO_MODEL - Names the JSON format file containing the
       PHYLO MODEL DATA; uses taxonomy.package.PHYLO_MODEL by default.
    """

    pkg_dir = args.package_name

    if args.clobber:
        utils.rmdir(pkg_dir)


    os.mkdir(pkg_dir)

    manifest = os.path.join(pkg_dir, package.MANIFEST_NAME)
    optdict = defaultdict(dict)

    # Add fields and values for the metadata section.
    optdict['metadata']['create_date'] = time.strftime('%Y-%m-%d %H:%M:%S')
    # locus is required
    optdict['metadata']['locus'] = args.locus
    # format_version is a module-level constant
    optdict['metadata']['format_version'] = package.FORMAT_VERSION
    if args.description:
        optdict['metadata']['description'] = args.description
    if args.author:
        optdict['metadata']['author'] = args.author
    if args.package_version:
        optdict['metadata']['package_version'] = args.package_version

    # phylo_model is part of the package, but not a command-line
    # argument; write out the phylo model file in JSON format, but
    # only if tree_stats was specified as an argument.
    if args.tree_stats:
        with open(os.path.join(pkg_dir, package.PHYLO_MODEL), 'w') as phylo_model:
            json.dump(utils.parse_raxml, phylo_model)

        optdict['files']['phylo_model'] = package.PHYLO_MODEL
        optdict['md5']['phylo_model'] = \
            hashlib.md5(open(phylo_model_pth).read()).hexdigest()

    # copy all provided files into the package directory
    for fname in package.PACKAGE_CONTENTS['files']:
        if fname == 'phylo_model':
            continue

        pth = getattr(args, fname)
        if pth:
            shutil.copy(pth, pkg_dir)
            optdict['files'][fname] = os.path.split(pth)[1]
            optdict['md5'][fname] = hashlib.md5(open(pth).read()).hexdigest()

    package.write_config(fname=manifest, optdict=optdict)


