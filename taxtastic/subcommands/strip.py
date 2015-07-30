"""Remove rollback and rollforward information from a refpkg

Delete everything in the refpkg not relevant to the current state,
including all files no longer referred to, as well as all rollback and
rollforward information. The log is preserved, with a new entry
entered indicating that ``refpkg`` was stripped.

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


def action(args):
    """Strips non-current files and rollback information from a refpkg.

    *args* should be an argparse object with fields refpkg (giving the
    path to the refpkg to operate on).
    """
    log.info('loading reference package')

    refpkg.Refpkg(args.refpkg, create=False).strip()
