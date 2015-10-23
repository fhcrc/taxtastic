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
"""Decide if tax_ids are primary and valid
"""

import csv
import logging
import pandas
import sys

import sqlalchemy

from sqlalchemy.sql import select
from taxtastic.taxonomy import Taxonomy
from taxtastic import ncbi

log = logging.getLogger(__name__)


def build_parser(parser):
    parser.add_argument(
        'infile',
        help="""Input CSV file to process, minimally containing the
        fields 'tax_id'. Rows with missing tax_ids are left unchanged.""")
    parser.add_argument(
        'database_file', help="""Path to the taxonomy database""")
    parser.add_argument(
        '-o', '--out-file', default=sys.stdout,
        help="""Output file to write updates [default: stdout]""")


def fetch_tax_info(tax_id, taxonomy):
    c = taxonomy.nodes.c  # columns
    s = select([c.tax_id, c.is_valid, c.parent_id, c.rank])
    s = s.where(c.tax_id == tax_id)
    return s.execute().fetchone()


def species_is_classified(tax_id, taxonomy):
    """
    returns is_valid from ncbi taxonomy
    """
    res = fetch_tax_info(tax_id, taxonomy)
    if not res:
        return False

    tax_id, is_valid, parent_id, rank = res
    if rank == 'species':
        return is_valid
    elif tax_id == parent_id:
        return False
    else:
        # instead of recursion create a table of ranks by index. code
        # will be something liek ranks.index(rank) > species_index
        return species_is_classified(parent_id, taxonomy)


def action(args):
    rows = pandas.read_csv(args.infile, dtype='str')

    if 'tax_id' not in rows.columns:
        raise ValueError("No tax_id column")

    if 'taxid_classified' in rows.columns:
        rows = rows.drop('taxid_classified', axis=1)

    e = sqlalchemy.create_engine('sqlite:///{0}'.format(args.database_file))
    tax = Taxonomy(e, ncbi.ranks)

    tax_ids = rows[['tax_id']].drop_duplicates()
    tax_ids.loc[:, 'taxid_classified'] = None

    def is_classified(row):
        row['taxid_classified'] = species_is_classified(row['tax_id'], tax)
        return row

    tax_ids = tax_ids.apply(is_classified, axis=1)
    rows = rows.merge(tax_ids, on='tax_id', how='left')
    rows.to_csv(args.out_file, index=False, quoting=csv.QUOTE_NONNUMERIC)
