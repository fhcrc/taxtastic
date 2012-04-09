"""
Adds or updates files or metdata in a refpkg.

The update subcommand takes a refpkg to operate on, then a series of changes to make, expressed as key=file.  So to add a file ../otherdir/boris under the key 'meep' and abcd under the key 'hilda' in a refpkg 'my-refpkg', you would run

$ taxit update my-refpkg meep=../otherdir/boris hilda=abcd

If a file already exists under a given key, it is overwritten.

Passing taxit update the --metadata option makes it update the metadata instead of files.  For example, to set the author field to "Genghis Khan" and the version to 0.4.3, run

$ taxit update --metadata "author=Genghis Khan" version=0.4.3
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

from taxtastic import refpkg

log = logging.getLogger(__name__)

def build_parser(parser):
    parser.add_argument('refpkg', action='store', metavar='refpkg',
                        help='the reference package to operate on')
    parser.add_argument('changes', nargs='*', metavar='key=value',
                        help='keys to update, in key=some_file format')
    parser.add_argument('--metadata', action='store_const', const=True,
                        default=False, help='Update metadata instead of files')


def action(args):
    """Updates a Refpkg with new files.

    *args* should be an argparse object with fields refpkg (giving the
    path to the refpkg to operate on) and changes (a series of strings
    of the form 'key=file' giving the key to update in the refpkg and
    the file to store under that key)."
    """
    log.info('loading reference package')

    pairs = [p.split('=',1) for p in args.changes]
    if args.metadata:
        rp = refpkg.Refpkg(args.refpkg, create=False)
        rp.start_transaction()
        for (key,value) in pairs:
            rp.update_metadata(key, value)
        rp.commit_transaction('Updated metadata: ' + \
                                  ', '.join(['%s=%s' % (a,b)
                                             for a,b in pairs]))
    else:
        for (key,filename) in pairs:
            if not(os.path.exists(filename)):
                print "No such file: %s" % filename
                exit(1)

        rp = refpkg.Refpkg(args.refpkg, create=False)
        rp.start_transaction()
        for (key,filename) in pairs:
            rp.update_file(key, os.path.abspath(filename))
        rp.commit_transaction('Updates files: ' + \
                                  ', '.join(['%s=%s' % (a,b)
                                             for a,b in pairs]))
    return 0
