"""
Rollforward a rolled back command on a refpkg.

$ taxit rollforward my-refpkg

You can also specify -n # to specify the number of operations to roll
forward (defaults to 1), as in

$ taxit rollforward -n 3 my-refpkg
"""

import logging
import os.path
import shutil
import sys

from taxtastic import refpkg

log = logging.getLogger(__name__)

def build_parser(parser):
    parser.add_argument('refpkg', action='store', metavar='refpkg',
                        help='the reference package to operate on')
    parser.add_argument('-n', action='store', metavar='int',
                        default=1, type=int, 
                        help='Number of operations to roll back')

def action(args):
    """Roll forward previously rolled back commands on a refpkg.

    *args* should be an argparse object with fields refpkg (giving the
    path to the refpkg to operate on) and optionall n (giving the
    number of operations to roll forward.
    """
    log.info('loading reference package')

    r = refpkg.Refpkg(args.refpkg)

    # First check if we can do n rollforwards
    q = r.contents
    for i in xrange(args.n):
        if q['rollforward'] == None:
            print >>sys.stderr, 'Cannot rollforward %d changes; refpkg only records %d rolled back changes.' % (args.n, i)
            return 1
        else:
            q = q['rollforward'][1]

    for i in xrange(args.n):
        r.rollforward()
    return 0

