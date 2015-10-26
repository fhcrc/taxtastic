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
"""Update obsolete tax_ids

Use in preparation for ``taxit taxtable``. Takes sequence info file as
passed to ``taxit create --seq-info``

"""

import csv
import logging
import pandas
import sys

import sqlalchemy

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
    parser.add_argument(
        '--unknowns',
        help="""csv file with single column 'seqname' identifying
        records with unknown taxids (works with --remove-unknown)""")
    parser.add_argument(
        '--name-column',
        help=('column with taxon name(s) to help '
              'find tax_ids. ex: organism name'))


def taxid_updater(taxonomy, halt=True):
    """
    Create a function to update obsolete taxonomies in
    """

    def update_taxid(tax_id, keyword=None):
        try:
            taxonomy._node(tax_id)
            # _node raises KeyError if the taxon couldn't be found in the
            # current taxonomy. If found, no update needed.
        except KeyError:
            new_tax_id = taxonomy._get_merged(tax_id)

            if new_tax_id != tax_id:
                msg = 'tax_id {} merged to {}'.format(tax_id, new_tax_id)
                log.info(msg)
                tax_id = new_tax_id
            elif keyword:
                try:
                    tax_id, _, _ = taxonomy.primary_from_name(keyword)
                except KeyError as err:
                    log.warn(err)
                    tax_id = None
            else:
                msg = "Unknown taxon {}".format(tax_id)
                log.warn(msg)
                tax_id = None
                if halt:
                    raise KeyError(msg)

        return tax_id

    return update_taxid


def action(args):
    rows = pandas.read_csv(args.infile, dtype='str')
    columns = rows.columns  # preserve column order

    if 'tax_id' not in columns:
        raise ValueError("No tax_id column")

    index = ['tax_id']

    if args.name_column:
        if args.name_column not in columns:
            raise ValueError("No " + args.name_column + " column")
        else:
            index.append(args.name_column)

    e = sqlalchemy.create_engine('sqlite:///{0}'.format(args.database_file))
    tax = Taxonomy(e, ncbi.ranks)

    tax_ids = rows[index].drop_duplicates()
    tax_ids.loc[:, 'new_tax_id'] = None

    updater = taxid_updater(tax, halt=not bool(args.unknowns))

    if args.name_column:
        def set_new_tax_id(row):
            row['new_tax_id'] = updater(
                row['tax_id'], keyword=row[args.name_column])
            return row
    else:
        def set_new_tax_id(row):
            row['new_tax_id'] = updater(row['tax_id'])
            return row

    tax_ids = tax_ids.apply(set_new_tax_id, axis=1)

    rows = rows.merge(tax_ids, on=index, how='left')
    rows = rows.drop('tax_id', axis=1)
    rows = rows.rename(columns={'new_tax_id': 'tax_id'})

    if args.unknowns:
        # unknown taxids are set to empty string in taxid_updater
        rows[rows['tax_id'].isnull()].to_csv(
            args.unknowns, index=False, columns=columns, quoting=csv.QUOTE_NONNUMERIC)

    rows.to_csv(args.out_file, index=False, columns=columns,
                quoting=csv.QUOTE_NONNUMERIC)
