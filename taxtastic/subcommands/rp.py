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
"""Resolve path; get the path to a file in the reference package

See online documentation for ``taxit create`` for a list of
permissible values for ``KEY``

For example, write the absolute path to the file containing the
phylogenetic tree in ``my.refpkg`` to stdout::

  taxit rp my.refpkg tree

Examine the contents of the seq_info file::

  less $(taxit rp my.refpkg seq_info)
"""
import logging
import sys

from taxtastic import refpkg

log = logging.getLogger(__name__)


def build_parser(parser):
    parser.add_argument('refpkg', action='store', metavar='refpkg',
                        help='the reference package to operate on')
    parser.add_argument('item', action='store', metavar='KEY',
                        help='show the path for file identified by KEY')


def action(args):
    rp = refpkg.Refpkg(args.refpkg, create=False)
    sys.stdout.write('%s\n' % rp.file_abspath(args.item))
    return 0
