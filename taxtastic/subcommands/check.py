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
"""Validate a reference package

Checks whether ``REFPKG`` is a valid input for ``pplacer``, that is,
does it have a FASTA file of the reference sequences; a Stockholm file
of their multiple alignment; a Newick formatted tree build from the
aligned sequences; and all the necessary auxiliary information.
"""

import taxtastic.refpkg


def build_parser(parser):
    parser.add_argument(
        'refpkg',
        action='store',
        metavar='REFPKG',
        help='Path to Refpkg to check')


def action(args):
    r = taxtastic.refpkg.Refpkg(args.refpkg, create=False)
    msg = r.is_ill_formed()
    if msg:
        print(msg)
        return 1
    else:
        return 0
