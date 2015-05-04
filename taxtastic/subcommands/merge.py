"""Identify merged taxids and provide replacements."""

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

import argparse
import csv
import logging

import re

from taxtastic import ncbi
from taxtastic.taxonomy import Taxonomy
from taxtastic.utils import getlines

from sqlalchemy import create_engine
import os.path
import sys

log = logging.getLogger(__name__)


def build_parser(parser):

    parser.add_argument(
        '-d', '--database-file',
        dest='database_file',
        metavar='FILE',
        required=True,
        help='Name of the sqlite database file')

    input_group = parser.add_argument_group(
        "Input options").add_mutually_exclusive_group()

    input_group.add_argument(
        '-t', '--tax-ids',
        dest='taxids',
        metavar='FILE-OR-LIST',
        help="""File containing a whitespace-delimited list of
        tax_ids (ie, separated by tabs, spaces, or newlines; lines
        beginning with "#" are ignored). This option can also be
        passed a comma-delited list of taxids on the command line.""")

    input_group.add_argument(
        '-i', '--seq-info', type=argparse.FileType('r'),
        help="""Read tax_ids from sequence info file, minimally containing the
        field "tax_id" """)

    output_group = parser.add_argument_group(
        "Output options").add_mutually_exclusive_group()

    output_group.add_argument(
        '-o', '--out-file',
        dest='out_file',
        type=argparse.FileType('w'),
        default=sys.stdout,
        metavar='FILE',
        help="""headerless csv file containing old_id,new_id""")


def action(args):
    engine = create_engine('sqlite:///%s' % args.database_file, echo=args.verbosity > 2)
    tax = Taxonomy(engine, ncbi.ranks)

    taxids = set()

    if args.taxids:
        if os.access(args.taxids, os.F_OK):
            for line in getlines(args.taxids):
                taxids.update(set(re.split(r'[\s,;]+', line)))
        else:
            taxids.update([x.strip() for x in re.split(r'[\s,;]+', args.taxids)])

    if args.seq_info:
        with args.seq_info:
            reader = csv.DictReader(args.seq_info)
            taxids.update(frozenset(i['tax_id'] for i in reader if i['tax_id']))

    writer = csv.writer(args.out_file)

    for t in taxids:
        try:
            tax._node(t)
        except KeyError:
            # Check for merged
            m = tax._get_merged(t)
            if m and m != t:
                writer.writerow([t, m])
            else:
                writer.writerow([t, None])

    engine.dispose()
    return 0
