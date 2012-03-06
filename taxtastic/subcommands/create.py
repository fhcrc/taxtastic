"""Creates a reference package"""
# This file is part of taxtastic.
#
#    taxtastic is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    taxtastic is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with taxtastic.  If not, see <http://www.gnu.org/licenses/>.

import logging
import shutil
import os
import time
import shutil
import hashlib
import re
import json
import sys

from taxtastic import refpkg

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
                        action="store", dest="seq_info", metavar="file",
                        help=('CSV format file describing the aligned reference '
                              'sequences, minimally containing the fields "seqname" '
                              'and "tax_id"'))
    parser.add_argument("-l", "--locus",
                        action="store", dest="locus", required=True,
                        help='The locus described by the reference package', metavar='LOCUS')
    parser.add_argument("-m", "--mask",
                        action="store", dest="mask",
                        help='Text file containing a mask', metavar='FILE')
    parser.add_argument("-p", "--profile",
                        action="store", dest="profile",
                        help='Alignment profile', metavar='FILE')
    parser.add_argument('-P', '--package-name', required=True,
                        action='store', dest='package_name',
                        metavar='PATH', help='Name of refpkg to create')
    parser.add_argument("-R", "--readme",
                        action="store", dest="readme",
                        help="README file describing the reference package",
                        metavar="FILE")
    parser.add_argument("-r", "--package-version",
                        action="store", dest="package_version",
                        help='Release version for the reference package', metavar='VERSION')
    parser.add_argument("-s", "--tree-stats",
                        action="store", dest="tree_stats",
                        help=('File containing tree statistics (for example '
                              'RAxML_info.whatever")'), metavar='FILE')
    parser.add_argument("-S", "--aln-sto",
                        action="store", dest="aln_sto",
                        help='Multiple alignment in Stockholm format', metavar='FILE')
    parser.add_argument("-t", "--tree-file",
                        action="store", dest="tree",
                        help='Phylogenetic tree in newick format',
                        metavar='FILE')
    parser.add_argument("-T", "--taxonomy",
                        action="store", dest="taxonomy",
                        help=('CSV format file defining the taxonomy. Fields include '
                              '"tax_id","parent_id","rank","tax_name" followed by a column '
                              'defining tax_id at each rank starting with root'),
                        metavar='FILE')


def action(args):
    if args.clobber and os.path.isdir(args.package_name):
        try:
            shutil.rmtree(args.package_name)
        except:
            print >>sys.stderr, "Failed: Could not delete %s" % args.package_name
            return 1
    elif args.clobber and os.path.exists(args.package_name):
        try:
            os.unlink(args.package_name)
        except:
            print >>sys.stderr, "Failed: Could not delete %s" % args.package_name
            return 1
    elif not args.clobber and os.path.exists(args.package_name):
        print >> sys.stderr, 'Failed: {0} exists.'.format(args.package_name)
        return 1

    r = refpkg.Refpkg(args.package_name)
    r.start_transaction()
    r.update_metadata('locus', args.locus) # Locus is required
    if args.description:
        r.update_metadata('description', args.description)
    if args.author:
        r.update_metadata('author', args.author)
    if args.package_version:
        r.update_metadata('package_version', args.package_version)
    if args.tree_stats:
        r.update_phylo_model(None, args.tree_stats)

    for file_name in ['aln_fasta', 'aln_sto', 'mask',
                      'profile', 'seq_info', 'taxonomy', 'tree', 'tree_stats',
                      'readme']:
        path = getattr(args, file_name)
        if path:
            r.update_file(file_name, path)
    r._log('Loaded initial files into empty refpkg')
    r.commit_transaction()
    r.strip()
    return 0
