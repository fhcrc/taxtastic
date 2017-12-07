"""Add or modify files or metadata in a refpkg

Update ``refpkg`` to set ``key`` to ``some value``.  If ``--metadata``
is specified, the update is done to the metadata.  Otherwise ``some
value`` is treated as the path to a file, and that file is updated in
``refpkg``.  An arbitrary of "key=value" pairs can be specified on the
command line.  If the same key is specified twice, the later
occurrence dominates.

All updates specified to an instance of this command are run as a
single operation, and will all be undone by a single rollback.

For example::

  taxit update my-refpkg meep=../otherdir/boris hilda=abcd

If a file already exists under a given key, it is overwritten.

The --metadata option causes a change to the metadata instead of
files.  For example, to set the author field to "Genghis Khan" and the
version to 0.4.3::

  taxit update --metadata "author=Genghis Khan" version=0.4.3

Other examples:

Set the author in my_refpkg::

  taxit update my_refpkg --metadata "author=Boris the mad baboon"

Set the author and version at once::

  taxit update my_refpkg --metadata "author=Bill" "package_version=1.7.2"

Insert a file into the refpkg::

  taxit update my_refpkg "aln_fasta=/path/to/a/file.fasta"

"""
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
import os.path
import warnings

from taxtastic import refpkg

log = logging.getLogger(__name__)


def build_parser(parser):
    parser.add_argument('refpkg', action='store', metavar='refpkg',
                        help='the reference package to operate on')
    parser.add_argument('changes', nargs='*', metavar='key=value',
                        help='keys to update, in key=some_file format')
    parser.add_argument('--metadata', action='store_const', const=True,
                        default=False, help='Update metadata instead of files')

    stats_group = parser.add_argument_group('Tree inference log file parsing '
                                            '(for updating `tree_stats`)')
    stats_group.add_argument("--stats-type", choices=('PhyML', 'FastTree', 'RAxML'),
                             help="""stats file type [default: attempt to guess from
                        file contents]""")
    stats_group.add_argument("--frequency-type", choices=('empirical', 'model'),
                             help="""Residue frequency type from the model. Required
                        for PhyML Amino Acid alignments.""")


def action(args):
    """Updates a Refpkg with new files.

    *args* should be an argparse object with fields refpkg (giving the
    path to the refpkg to operate on) and changes (a series of strings
    of the form 'key=file' giving the key to update in the refpkg and
    the file to store under that key)."
    """
    log.info('loading reference package')

    pairs = [p.split('=', 1) for p in args.changes]
    if args.metadata:
        rp = refpkg.Refpkg(args.refpkg, create=False)
        rp.start_transaction()
        for key, value in pairs:
            rp.update_metadata(key, value)
        rp.commit_transaction('Updated metadata: ' +
                              ', '.join(['%s=%s' % (a, b)
                                         for a, b in pairs]))
    else:
        for key, filename in pairs:
            if not(os.path.exists(filename)):
                print("No such file: %s" % filename)
                exit(1)

        rp = refpkg.Refpkg(args.refpkg, create=False)
        rp.start_transaction()
        for key, filename in pairs:
            if key == 'tree_stats':
                with warnings.catch_warnings():
                    warnings.simplefilter(
                        "ignore", refpkg.DerivedFileNotUpdatedWarning)
                    rp.update_file(key, os.path.abspath(filename))
                # Trigger model update
                log.info('Updating phylo_model to match tree_stats')
                rp.update_phylo_model(args.stats_type, filename,
                                      args.frequency_type)
            else:
                rp.update_file(key, os.path.abspath(filename))

        rp.commit_transaction('Updates files: ' +
                              ', '.join(['%s=%s' % (a, b)
                                         for a, b in pairs]))
    return 0
