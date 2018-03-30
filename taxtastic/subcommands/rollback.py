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
"""Undo an operation performed on a refpkg

Rollback ``N`` operations on ``refpkg`` (default to 1 operation if
``-n`` is omitted).  This is equivalent to calling the ``rollback()``
method of ``taxtastic.refpkg.Refpkg``.  If there are not at least
``N`` operations that can be rolled back, an error is returned and no
changes are made to the refpkg.
"""

import logging

from taxtastic import refpkg

log = logging.getLogger(__name__)


def build_parser(parser):
    parser.add_argument('refpkg', action='store', metavar='refpkg',
                        help='the reference package to operate on')
    parser.add_argument('-n', action='store', metavar='int',
                        default=1, type=int,
                        help='Number of operations to roll back')


def action(args):
    """Roll back commands on a refpkg.

    *args* should be an argparse object with fields refpkg (giving the
    path to the refpkg to operate on) and n (giving the number of
    operations to roll back).
    """
    log.info('loading reference package')

    r = refpkg.Refpkg(args.refpkg, create=False)

    # First check if we can do n rollbacks
    q = r.contents
    for i in range(args.n):
        if q['rollback'] is None:
            log.error('Cannot rollback {} changes; '
                      'refpkg only records {} changes.'.format(args.n, i))
            return 1
        else:
            q = q['rollback']

    for i in range(args.n):
        r.rollback()

    return 0
