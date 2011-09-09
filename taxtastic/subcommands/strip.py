"""
Removes all rollback and rollforward information and files not attached to the current state from a refpkg.

$ taxit strip my-refpkg
"""

import logging
import os.path
import shutil

from taxtastic import refpkg

log = logging.getLogger(__name__)

def build_parser(parser):
    parser.add_argument('refpkg', action='store', metavar='refpkg',
                        help='the reference package to operate on')


def action(args):
    """Strips non-current files and rollback information from a refpkg.

    *args* should be an argparse object with fields refpkg (giving the
    path to the refpkg to operate on).
    """
    log.info('loading reference package')

    refpkg.Refpkg(args.refpkg).strip()
