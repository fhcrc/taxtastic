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
"""Create a tabular representation of taxonomic lineages

Write a CSV file containing the minimal subset of the taxonomy in
``database_file`` representing all of the lineages specified by the
provided tax_ids. Duplicate tax_ids are ignored.

By default the CSV is written to ``stdout``, unless a file is
specified with ``-o/--outfile``.

"""
import argparse
import csv
import logging
import re
import sqlalchemy
import sys
from itertools import groupby

from taxtastic.taxonomy import Taxonomy
from taxtastic.utils import add_database_args

log = logging.getLogger(__name__)


def replace_no_rank(ranks):
    ranks = list(ranks)
    while True:
        try:
            idx = ranks.index('no_rank')
            ranks[idx] = ranks[idx - 1] + '_'
        except ValueError:
            return ranks


def as_taxtable_rows(rows, seen=None):

    # allows termination when we've seen a tax_id before - shaves just
    # a few seconds off the total time.
    seen = seen or {}

    __, tids, pids, ranks, names = [list(tup) for tup in zip(*rows)]
    ranks = replace_no_rank(ranks)
    ranks_out = ranks[:]

    tax_rows = []
    while tids and tids[-1] not in seen:
        d = dict(zip(ranks, tids))
        d['tax_id'] = tids.pop(-1)
        d['parent_id'] = pids.pop(-1)
        d['tax_name'] = names.pop(-1)
        d['rank'] = ranks.pop(-1)
        tax_rows.append((d['tax_id'], d))

    return ranks_out, tax_rows


def order_ranks(ref_ranks):
    def _inner(rank):
        trailing_ = re.findall(r'_+$', rank)
        if trailing_:
            return (ref_ranks.index(rank.rstrip('_')), len(trailing_[0]))
        else:
            return (ref_ranks.index(rank), 0)

    return _inner


def getitems(*items):
    def _getter(d):
        return [d.get(item) for item in items]

    return _getter


def build_parser(parser):
    parser = add_database_args(parser)

    # TODO: do we need the commented-out options?
    # node_parser = parser.add_argument_group(title='node options')

    # node_parser.add_argument(
    #     '--valid',
    #     action='store_true',
    #     help='Include only valid nodes.')
    # node_parser.add_argument(
    #     '--ranked',
    #     choices=['columns', 'rows'],
    #     help='Include only ranked columns or ranked rows and columns.')

    input_group = parser.add_argument_group('input options')

    # input_group.add_argument(
    #     '--taxtable',
    #     metavar='CSV',
    #     help='build from a previous built taxtable')

    # input_group.add_argument(
    #     '--clade-ids',
    #     help=('return top-down tax_id clades'))

    input_group.add_argument(
        '-t', '--tax-ids', nargs='+',
        help='one or more tax_ids')

    input_group.add_argument(
        '-f', '--tax-id-file', metavar='FILE', type=argparse.FileType('r'),
        help=('File containing a whitespace-delimited list of '
              'tax_ids (ie, separated by tabs, spaces, or newlines.'))

    input_group.add_argument(
        '-i', '--seq-info',
        type=argparse.FileType('r'),
        help=('Read tax_ids from sequence info file, minimally '
              'containing a columnm named "tax_id"'))

    output_group = parser.add_argument_group(
        "Output options").add_mutually_exclusive_group()

    output_group.add_argument(
        '-o', '--outfile',
        type=argparse.FileType('w'),
        default=sys.stdout,
        metavar='FILE',
        help=('Output file containing lineages for the specified taxa '
              'in csv format; writes to stdout if unspecified'))


def action(args):
    log.info('reading tax_ids')
    if args.tax_ids:
        tax_ids = set(args.tax_ids)
    elif args.tax_id_file:
        tax_ids = set(args.tax_id_file.read().split())
    elif args.seq_info:
        tax_ids = {row['tax_id'] for row in csv.DictReader(args.seq_info)}
    else:
        sys.exit('Error: no tax_ids were specified')

    engine = sqlalchemy.create_engine(args.url, echo=args.verbosity > 3)
    tax = Taxonomy(engine, schema=args.schema)

    rows = tax._get_lineage_table(tax_ids)
    log.info('grouping lineages')
    all_ranks = set()
    taxtable = {}
    for tax_id, grp in groupby(rows, lambda row: row[0]):
        ranks, tax_rows = as_taxtable_rows(grp, seen=taxtable)
        taxtable.update(dict(tax_rows))
        all_ranks |= set(ranks)

    sorted_ranks = sorted(all_ranks, key=order_ranks(tax.ranks[::-1]))
    fieldnames = ['tax_id', 'parent_id', 'rank', 'tax_name'] + sorted_ranks
    output = taxtable.values()
    log.info('sorting lineages')
    output = sorted(output, key=getitems(*sorted_ranks))
    # Now do a sanity check to be sure every row has a parent_id. If the parent node is empty (ie root) make it self-referential 
    # There are fancier ways to do this but a for loop will suffice. I'll use indices to replace in-situ and preserve memory
    for i in xrange(len(output)):
        if output[i]['parent_id'] == None or output[i]['parent_id'] == "":
            output[i]['parent_id'] = output[i]['tax_id']

    log.info('writing taxtable')
    writer = csv.DictWriter(
        args.outfile, fieldnames=fieldnames, quoting=csv.QUOTE_ALL)
    writer.writeheader()
    writer.writerows(output)
