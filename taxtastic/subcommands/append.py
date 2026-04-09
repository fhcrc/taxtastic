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
"""Append taxonomy columns to a CSV file.

Given a CSV file with a tax_id column, look up the taxonomic lineage
for each tax_id and append the requested rank columns (e.g. genus,
species) to each row.  Use _name to for rank tax_name (e.g. genus_name,
species_name)

Example::

    taxit append --columns genus,species,species_name seq_info.csv ncbi.db

"""
import csv
import logging
import sys
from collections import defaultdict

import sqlalchemy

from taxtastic.taxonomy import Taxonomy
from taxtastic.utils import Opener, add_database_args

log = logging.getLogger(__name__)


def build_parser(parser):
    parser.add_argument(
        'infile',
        type=Opener('r'),
        help='input CSV file with a tax_id column')
    parser = add_database_args(parser)
    parser.add_argument(
        '-c', '--columns',
        type=lambda x: [c.strip() for c in x.split(',')],
        required=True,
        metavar='COLUMNS',
        help='comma-separated list of taxonomy ranks and tax_names '
             'to append (e.g. genus,genus_name,species,species_name).')
    parser.add_argument(
        '--tax-id-column',
        default='tax_id',
        metavar='COLUMN',
        help='name of the tax_id column in the input CSV [%(default)s]')
    parser.add_argument(
        '-o', '--outfile',
        type=Opener('w'),
        default=sys.stdout,
        metavar='FILE',
        help='output CSV file [stdout]')


def action(args):
    engine = sqlalchemy.create_engine(args.url, echo=args.verbosity > 3)
    tax = Taxonomy(engine, schema=args.schema)

    reader = csv.DictReader(args.infile)
    rows = list(reader)

    tax_ids = {
        row[args.tax_id_column]
        for row in rows
        if args.tax_id_column in row
    }

    log.info('fetching lineages for {} tax_ids'.format(len(tax_ids)))
    lineages = defaultdict(dict)
    try:
        lineage_rows = tax._get_lineage_table(tax_ids)
        for tid, _tax_id, _parent_id, rank, tax_name in lineage_rows:
            lineages[tid][rank] = (_tax_id, tax_name)
    except ValueError:
        log.warning('no tax_ids were found in the database')

    missing = tax_ids - set(lineages)
    if missing:
        log.warning(
            '{} tax_ids not found in database: {}'.format(
                len(missing), sorted(missing)))

    existing = set(reader.fieldnames or [])
    new_cols = [c for c in args.columns if c not in existing]
    fieldnames = (reader.fieldnames or []) + new_cols

    writer = csv.DictWriter(args.outfile, fieldnames=fieldnames)
    writer.writeheader()
    for row in rows:
        tid = row.get(args.tax_id_column, '')
        lineage = lineages.get(tid, {})
        for col in args.columns:
            rank = col.rstrip('_name')
            _tax_id, tax_name = lineage.get(rank, ('', ''))
            if col.endswith('_name'):
                row[col] = tax_name
            else:
                row[col] = _tax_id
        writer.writerow(row)

    engine.dispose()
