"""
Rollforward a rolled back command on a refpkg.

$ taxit rollforward my-refpkg

You can also specify -n # to specify the number of operations to roll
forward (defaults to 1), as in

$ taxit rollforward -n 3 my-refpkg
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

    r = refpkg.Refpkg(args.refpkg, create=False)

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

