"""Updates the contents of a reference package."""

import logging
import os.path
import shutil

from taxtastic import refpkg

log = logging.getLogger(__name__)

def build_parser(parser):
    parser.add_argument('refpkg', nargs=1,
        help='the reference package to operate on')
    parser.add_argument('changes', nargs='*', metavar='key=value',
        help='keys to update, in key=some_file format')

def action(args):
    log.info('loading reference package')
    rp = refpkg.Refpkg(args.refpkg[0])

    for pair in args.changes:
        key, _, value = pair.partition('=')
        outfile_name = os.path.basename(value)
        with open(value, 'rb') as infile, \
                rp.file_resource(outfile_name, 'wb') as outfile:
            shutil.copyfileobj(infile, outfile)
        rp.contents['files'][key] = outfile_name
        rp.rehash(key)

    rp.save()
