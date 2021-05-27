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
"""Create a reference package

Create a new refpkg at the location specified by the argument to
``-P`` with locus name ``-l``.  All other fields are used to specify
initial metadata and files to add to the refpkg.  If there is already
a refpkg at ``refpkg``, this command will fail unless you specify
``-c`` or ``--clobber``."""
import logging
import shutil
import os

from taxtastic import refpkg, utils

log = logging.getLogger(__name__)


class ConfigError(Exception):
    pass


def build_parser(parser):

    parser.add_argument("-c", "--clobber",
                        action="store_true", dest="clobber", default=False,
                        help='Delete an existing reference package.')

    required = parser.add_argument_group('Required arguments')
    required.add_argument('-P', '--package-name', required=True,
                          action='store', dest='package_name',
                          metavar='PATH', help='Name of refpkg to create')
    required.add_argument("-l", "--locus",
                          action="store", dest="locus", required=True,
                          help='The locus described by the reference package',
                          metavar='LOCUS')

    metadata = parser.add_argument_group('Package Metadata')

    metadata.add_argument("-a", "--author",
                          action="store", dest="author",
                          help='Person who created the reference package', metavar='NAME')
    metadata.add_argument("-d", "--description",
                          action="store", dest="description",
                          help='An arbitrary description field', metavar='TEXT')
    metadata.add_argument("-r", "--package-version",
                          action="store", dest="package_version",
                          help='Release version for the reference package',
                          metavar='VERSION')

    infiles = parser.add_argument_group('Input files')

    infiles.add_argument("-f", "--aln-fasta",
                         action="store", dest="aln_fasta",
                         help='Multiple alignment in fasta format', metavar='FILE')
    infiles.add_argument("-i", "--seq-info",
                         action="store", dest="seq_info", metavar="file",
                         help=('CSV format file describing the aligned reference '
                               'sequences, minimally containing the fields "seqname" '
                               'and "tax_id"'))
    infiles.add_argument("-m", "--mask",
                         action="store", dest="mask",
                         help='Text file containing a mask', metavar='FILE')
    infiles.add_argument("-M", "--model-file", action="store", dest="model",
                         help=('File containing model information '
                               'usually the .bestModel file'), metavar='FILE')
    infiles.add_argument("-p", "--profile",
                         action="store", dest="profile",
                         help='Alignment profile', metavar='FILE')
    infiles.add_argument("-R", "--readme",
                         action="store", dest="readme",
                         help="README file describing the reference package",
                         metavar="FILE")
    infiles.add_argument("-s", "--tree-stats",
                         action="store", dest="tree_stats",
                         help=('File containing tree statistics (for example '
                               'RAxML_info.whatever")'), metavar='FILE')
    infiles.add_argument("-S", "--aln-sto",
                         action="store", dest="aln_sto",
                         help='Multiple alignment in Stockholm format', metavar='FILE')
    infiles.add_argument("-t", "--tree-file",
                         action="store", dest="tree",
                         help='Phylogenetic tree in newick format',
                         metavar='FILE')
    infiles.add_argument(
        "-T", "--taxonomy",
        action="store", dest="taxonomy",
        help=('CSV format file defining the taxonomy. Fields include '
              '"tax_id","parent_id","rank","tax_name" followed by a column '
              'defining tax_id at each rank starting with root'),
        metavar='FILE')

    tree_info = parser.add_argument_group('Tree information')
    tree_info.add_argument("--stats-type", choices=('PhyML', 'FastTree', 'RAxML'),
                           help="""stats file type [default: attempt to guess from
                           file contents]""")
    tree_info.add_argument(
        "--frequency-type", choices=('empirical', 'model'),
        help="""Residue frequency type from the model. Required for
        var in collection: PhyML Amino Acid alignments.""")

    root_grp = parser.add_argument_group('Taxonomic Rerooting')
    root_grp.add_argument(
        '--no-reroot', action='store_false', dest='reroot',
        default=True, help="""Do not reroot the reference package
        using `rppr reroot`. [default: reroot if `rppr` is available
        and a taxonomy file is specified]""")
    root_grp.add_argument('--rppr', default='rppr', help="""Name of the rppr
            executable. [default: %(default)s]""")


def action(args):
    if args.clobber and os.path.isdir(args.package_name):
        try:
            shutil.rmtree(args.package_name)
        except Exception:
            log.error("Failed: Could not delete %s" % args.package_name)
            return 1
    elif args.clobber and os.path.exists(args.package_name):
        try:
            os.unlink(args.package_name)
        except Exception:
            log.error("Failed: Could not delete %s" % args.package_name)
            return 1
    elif not args.clobber and os.path.exists(args.package_name):
        log.error('Failed: {0} exists.'.format(args.package_name))
        return 1

    r = refpkg.Refpkg(args.package_name, create=True)
    r.start_transaction()
    r.update_metadata('locus', args.locus)  # Locus is required
    if args.description:
        r.update_metadata('description', args.description)
    if args.author:
        r.update_metadata('author', args.author)
    if args.package_version:
        r.update_metadata('package_version', args.package_version)
    if args.tree_stats:
        r.update_phylo_model(args.stats_type, args.tree_stats,
                             frequency_type=args.frequency_type)

    for file_name in ['aln_fasta', 'aln_sto', 'mask',
                      'profile', 'seq_info', 'taxonomy', 'tree', 'tree_stats',
                      'readme', 'model']:
        path = getattr(args, file_name)
        if path:
            r.update_file(file_name, path)
    r._log('Loaded initial files into empty refpkg')
    r.commit_transaction()
    r.strip()

    reroot_prereqs = args.reroot and args.taxonomy and args.seq_info and args.tree
    if utils.has_rppr(args.rppr) and reroot_prereqs:
        r.start_transaction()
        logging.info('%s found. Rerooting.', args.rppr)
        r.reroot(rppr=args.rppr)
        r._log('Rerooted')
        r.commit_transaction()
    elif reroot_prereqs:
        log.warn('"%s" not found. Skipping rerooting', args.rppr)

    return 0
