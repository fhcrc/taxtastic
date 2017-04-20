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
"""Extracts tax ids of all lonely nodes in a taxtable

Find nodes in ``target`` (which can be a CSV file extracted by ``taxit
taxtable`` or a RefPkg containing such a file) which are lonely; that
is, whose parents have only one child. Print them, one per line, to
``stdout`` or to the file specified by the ``-o`` option.
"""
import argparse
import logging
import os
import sys
import csv

from taxtastic import lonely, refpkg

log = logging.getLogger(__name__)


class ConfigError(Exception):
    pass


def comma_separated_set(s):
    return frozenset(i.strip() for i in s.split(','))


def build_parser(parser):
    parser.add_argument("target",
                        metavar="taxtable_or_refpkg",
                        action="store",
                        help='A taxtable or a refpkg containing a taxtable')
    parser.add_argument('-o', '--out',
                        type=argparse.FileType('w'),
                        default=sys.stdout,
                        help="""Write output to given file [default: stdout]""")
    parser.add_argument('-r', '--ranks', help="""Comma separated list of ranks
            to consider [default: all ranks]""", type=comma_separated_set)


def action(args):
    if not(os.path.exists(args.target)):
        logging.error("Failed: no such target %s", args.target)
        return 1
    elif os.path.isdir(args.target):
        logging.info("Target is a refpkg. Working on taxonomy within it.")
        r = refpkg.Refpkg(args.target, create=False)
        path = r.file_abspath('taxonomy')
    else:
        logging.info("Target is a CSV file")
        path = args.target

    logging.info("Loading taxonomy from file.")
    with open(path, 'rU') as h:
        tree = lonely.taxtable_to_tree(h)
    result = tree.lonelynodes()
    if args.ranks:
        result = (n for n in result if n.rank in args.ranks)

    writer = csv.writer(args.out)
    writer.writerow(['tax_name', 'tax_id', 'rank'])
    writer.writerows(sorted((n.tax_name, n.key, n.rank) for n in result))
