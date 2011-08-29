"""
Adds or updates files or metdata in a refpkg.

The update subcommand takes a refpkg to operate on, then a series of changes to make, expressed as key=file.  So to add a file ../otherdir/boris under the key 'meep' and abcd under the key 'hilda' in a refpkg 'my-refpkg', you would run

$ taxit update my-refpkg meep=../otherdir/boris hilda=abcd

If a file already exists under a given key, it is overwritten.

Passing taxit update the --metadata option makes it update the metadata instead of files.  For example, to set the author field to "Genghis Khan" and the version to 0.4.3, run

$ taxit update --metadata "author=Genghis Khan" version=0.4.3
"""

import logging
import os.path
import shutil

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
        rp = refpkg.Refpkg(args.refpkg)
        for (key,value) in pairs:
            rp.update_metadata(key, value)
    else:
        for (key,filename) in pairs:
            if not(os.path.exists(filename)):
                print "No such file: %s" % filename
                exit(1)

        rp = refpkg.Refpkg(args.refpkg)
        for (key,filename) in pairs:
            rp.update_file(key, os.path.abspath(filename))
