"""Taxonomically reroots a reference package

Calls ``rppr reroot`` to generate a rerooted tree from the tree in
``refpkg`` and writes it back to the refpkg.  The refpkg ``refpkg``
must contain the necessary inputs for ``pplacer`` for this to work.

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

from taxtastic import refpkg

log = logging.getLogger(__name__)


def build_parser(parser):
    parser.add_argument('refpkg', action='store', metavar='refpkg',
                        help='the reference package to operate on')
    parser.add_argument('--rppr', action='store', default=None,
                        help="specify the rppr binary to call to perform the rerooting")
    parser.add_argument('-p', '--pretend',
                        action='store_true', default=False,
                        help="don't save the rerooted tree; just attempt the rerooting.")


def action(args):
    r = refpkg.Refpkg(args.refpkg, create=False)
    r.reroot(rppr=args.rppr, pretend=args.pretend)
