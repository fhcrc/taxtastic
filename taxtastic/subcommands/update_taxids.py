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

FIXME: only include --name-column column on unknown tax_ids.  As is we are
losing a lot of efficiency by dedupicating on both tax_id,organism columns
for all rows.

TODO: implement --append-column RANK.
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
        nargs='?',
        default=sys.stdin,
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
        '--taxid-classified', action='store_true',
        help='add column True/False column if the tax_id is primary and valid')
    parser.add_argument(
        '--name-column',
        help=('column with taxon name(s) to help '
              'find tax_ids. ex: organism name'))
    parser.add_argument(
        '--append-lineage',
        help=('rank to append to seq_info'))


def update_taxid(tax_id, taxonomy, halt=True):
    try:
        taxonomy._node(tax_id)
        # _node raises ValueError if the taxon couldn't be found in the
        # current taxonomy. If found, no update needed.
    except ValueError as err:
        new_tax_id = taxonomy._get_merged(tax_id)

        if new_tax_id != tax_id:
            log.info('tax_id {} merged to {}'.format(tax_id, new_tax_id))
            tax_id = new_tax_id
        else:
            if halt:
                raise err
            else:
                log.warn(err)
                tax_id = None

    return tax_id


def update_by_name(name, taxonomy, halt=True):
    try:
        tax_id, _, _ = taxonomy.primary_from_name(name)
    except ValueError as err:
        if halt:
            raise err
        else:
            log.warn(err)
            tax_id = None
    return tax_id


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


def fetch_tax_info(tax_id, taxonomy):
    c = taxonomy.nodes.c  # columns
    s = select([c.tax_id, c.is_valid, c.parent_id, c.rank])
    s = s.where(c.tax_id == tax_id)
    return s.execute().fetchone()


def action(args):
    rows = pandas.read_csv(args.infile, dtype='str')
    columns = rows.columns.tolist()  # preserve column order

    if 'tax_id' not in columns:
        raise ValueError("No tax_id column")

    if args.name_column:
        if args.name_column not in columns:
            raise ValueError("No " + args.name_column + " column")

    e = sqlalchemy.create_engine('sqlite:///{0}'.format(args.database_file))
    tax = Taxonomy(e, ncbi.ranks)

    row_total = len(rows)

    def print_status(msg, row_num):
        sys.stderr.write(msg + ' {:.2%}\r'.format(float(row_num) / row_total))

    tax_ids = rows['tax_id'].drop_duplicates()
    halt = not (bool(args.unknowns) or bool(args.name_column))

    def set_new_tax_id(row):
        print_status('updating taxids: ', row.name)
        row['new_tax_id'] = update_taxid(row['tax_id'], tax, halt=halt)
        return row

    tax_ids = tax_ids.apply(set_new_tax_id, axis=1)

    if args.name_column:
        rows = rows.merge(tax_ids, on='tax_id', how='left')
        no_ids = rows['new_tax_id'].isnull()
        rows_no_ids = rows[no_ids].drop('new_tax_id', axis=1)

        name_ids = rows_no_ids[args.name_column].drop_duplicates()

        halt = not bool(args.unknowns)

        def set_name_taxid(row):
            row['new_tax_id'] = update_by_name(row['name'], tax, halt=halt)
            return row

        log.info('searching unknown tax_ids by ' + args.name_column)
        name_ids = name_ids.apply(set_name_taxid, axis=1)
        rows_no_ids = rows_no_ids.merge(
            name_ids, on=args.name_column, how='left')

        rows = pandas.concat(rows[~no_ids], rows_no_ids)
        tax_ids = rows['tax_id'].drop_duplicates()

    if args.taxid_classified:
        if 'taxid_classified' in columns:
            rows = rows.drop('taxid_classified', axis=1)
        else:
            columns.append('taxid_classified')

        tax_ids.loc[:, 'taxid_classified'] = None

        def is_classified(row):
            print_status('determining taxon name validity: ', row.name)
            row['taxid_classified'] = species_is_classified(
                row['tax_id'], tax)
            return row

        tax_ids = tax_ids.apply(is_classified, axis=1)

    if args.append_lineage:
        if args.append_lineage in columns:
            rows = rows.drop(args.append_lineage, axis=1)
        else:
            columns.append(args.append_lineage)

        tax_ids.loc[:, args.append_lineage] = None

        msg = 'appending lineage {}: '.format(args.append_lineage)

        def add_rank_column(row):
            print_status(msg, row.name)
            lineage = tax.lineage(row['tax_id'])
            row[args.append_lineage] = lineage.get(args.append_lineage, None)
            return row

        tax_ids = tax_ids.apply(add_rank_column, axis=1)

    rows = rows.merge(tax_ids, on='tax_id', how='left')

    if args.unknowns:
        # unknown taxids are set to empty string in taxid_updater
        rows[rows['tax_id'].isnull()].to_csv(
            args.unknowns, index=False,
            columns=columns, quoting=csv.QUOTE_NONNUMERIC)

    rows.to_csv(args.out_file, index=False, columns=columns,
                quoting=csv.QUOTE_NONNUMERIC)
