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
"""Add nodes to an existing taxtable csv"""
import argparse
import csv
import logging
import sys

from taxtastic import taxtable

log = logging.getLogger(__name__)


def build_parser(parser):
    parser.add_argument(
        "taxtable",
        metavar='CSV',
        type=argparse.FileType('r'),
        help="""A taxtable to augment""")
    parser.add_argument(
        'extra_nodes_csv',
        metavar='CSV',
        type=argparse.FileType('r'),
        help="""A CSV
        file containing nodes to add to taxtable. Must contain columns
        'tax_id', 'tax_name', 'rank', and 'parent_id'. Each record
        must have a parent_id already in the taxtable, or defined on
        an earlier row.""")
    parser.add_argument(
        '-o', '--out',
        type=argparse.FileType('w'),
        default=sys.stdout,
        metavar='CSV',
        help="""Destination for output taxtable [default: stdout]""")


def add_rank(tax, parent_node, rank):
    parent_idx = tax.ranks.index(parent_node.rank)
    assert parent_idx >= 0, "Parent rank must exist"

    # Insert after parent rank.
    # This isn't perfect, since there's no guarantee that there aren't
    # intervening ranks, but we don't have enough information to sort that.
    tax.ranks.insert(parent_idx + 1, rank)


def action(args):
    with args.taxtable as fp:
        tax = taxtable.read(fp)

    with args.extra_nodes_csv:
        reader = csv.DictReader(args.extra_nodes_csv)
        missing_fields = frozenset(
            ['tax_id', 'tax_name', 'rank', 'parent_id']) - frozenset(reader.fieldnames)
        if missing_fields:
            raise IOError("Missing expected fields: {0}".format(
                ','.join(missing_fields)))
        for row in reader:
            if row['tax_id'] in tax.index:
                logging.warn("tax_id %s already represented in taxtable. [row %d]",
                             row['tax_id'], reader.line_num)
                continue

            parent_id = row['parent_id']
            rank = row['rank']
            try:
                parent_node = tax.get_node(parent_id)
            except ValueError:
                raise ValueError(
                    "Parent {parent_id} of {tax_id}[{tax_name}] not found.".format(**row))
            if rank not in tax.ranks:
                add_rank(tax, parent_node, rank)
            node = taxtable.TaxNode(
                tax_id=row['tax_id'], name=row['tax_name'], rank=rank)
            parent_node.add_child(node)
            logging.info(
                "Added %s %s[%s] below %s %s[%s]",
                node.rank, node.tax_id, node.name,
                parent_node.rank, parent_node.tax_id, parent_node.name)

    tax.write_taxtable(args.out)
