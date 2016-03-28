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

TODO: refactor to simplify the action function.
"""

import csv
import logging
import pandas
import sqlalchemy
import sys

from sqlalchemy.sql import select

from taxtastic.taxonomy import Taxonomy
from taxtastic import ncbi, utils

log = logging.getLogger(__name__)


def build_parser(parser):
    parser.add_argument(
        'infile',
        nargs='?',
        default=sys.stdin,
        help=('Input CSV file to process, minimally '
              'containing the fields `tax_id`. Rows with '
              'missing tax_ids are left unchanged.'))
    parser.add_argument(
        'database_file', help="""Path to the taxonomy database""")
    parser.add_argument(
        '-o', '--out-file',
        default=sys.stdout,
        help='Output file to write updates [default: stdout]')
    parser.add_argument(
        '--taxid-column',
        default='tax_id',
        help='name of tax_id column to update [%(default)s]')
    parser.add_argument(
        '--unknowns',
        help=('unchanged taxonomy output of records with unknown taxids'))
    parser.add_argument(
        '--ignore-unknowns',
        action='store_true',
        help='allow unknown tax_ids in final output')
    parser.add_argument(
        '--taxid-classified',
        action='store_true',
        help=('add column True/False column if '
              'the tax_id is primary and valid'))
    parser.add_argument(
        '--name-column',
        help=('column with taxon name(s) to help '
              'find tax_ids. ex: organism name'))
    parser.add_argument(
        '--append-lineage',
        help=('rank to append to seq_info'))


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

    if args.taxid_column not in columns:
        raise ValueError("No column " + args.taxid_column)

    if args.name_column:
        if args.name_column not in columns:
            msg = '"No "' + args.name_column + '" column'
            raise ValueError(msg)

    con = 'sqlite:///{0}'.format(args.database_file)
    e = sqlalchemy.create_engine(con)
    tax = Taxonomy(e, ncbi.ranks)

    merged = pandas.read_sql_table('merged', con, index_col='old_tax_id')
    log.info('updating tax_ids')
    rows = rows.join(merged, on=args.taxid_column)

    # overwrite tax_ids where there is a new_tax_id
    inew_tax_ids = ~rows['new_tax_id'].isnull()
    rows.loc[inew_tax_ids, args.taxid_column] = \
        rows[inew_tax_ids]['new_tax_id']
    rows = rows.drop('new_tax_id', axis=1)

    log.info('loading names table')
    names = pandas.read_sql_table(
        'names', con, columns=['tax_id', 'tax_name', 'is_primary'])

    if args.name_column:
        """
        use the args.name_column to do a string comparison with
        names.tax_name column to find a suitable tax_id
        """
        unknowns = rows[~rows[args.taxid_column].isin(names['tax_id'])]

        if not unknowns.empty:
            """
            Take any tax_id associated with a string match
            to tax_name prioritizing is_primary=True
            """
            unknowns = unknowns.drop(args.taxid_column, axis=1)
            names = names.sort_values('is_primary', ascending=False)
            names = names.drop_duplicates(subset='tax_name', keep='first')
            names = names.set_index('tax_name')
            found = unknowns.join(names, on=args.name_column, how='inner')
            rows.loc[found.index, args.taxid_column] = found['tax_id']

    if not args.ignore_unknowns:
        unknowns = rows[~rows[args.taxid_column].isin(names['tax_id'])]
        if args.unknowns:
            """
            Output unknown tax_ids
            """
            unknowns.to_csv(
                args.unknowns,
                index=False,
                columns=columns,
                quoting=csv.QUOTE_NONNUMERIC)
        elif not unknowns.empty:
            raise ValueError('Unknown or missing tax_ids present')

        rows = rows[~rows.index.isin(unknowns.index)]

    if args.taxid_classified:
        """
        """
        if 'taxid_classified' in columns:
            rows = rows.drop('taxid_classified', axis=1)
        else:
            columns.append('taxid_classified')

        def is_classified(row):
            row['taxid_classified'] = species_is_classified(
                row[args.taxid_column], tax)
            return row

        msg = 'validating tax_ids:'
        rows = utils.apply_df_status(is_classified, rows, msg)

    if args.append_lineage:
        """
        Append a column from the taxonomy to seq_info
        """
        if args.append_lineage in columns:
            rows = rows.drop(args.append_lineage, axis=1)
        else:
            columns.append(args.append_lineage)

        def add_rank_column(row):
            try:
                lineage = tax.lineage(row[args.taxid_column])
            except ValueError as e:
                log.warn(e)
                lineage = {}
            row[args.append_lineage] = lineage.get(args.append_lineage, None)
            return row

        msg = 'appending {} column'.format(args.append_lineage)
        rows = utils.apply_df_status(add_rank_column, rows, msg)

    # output seq_info with new tax_ids
    rows.to_csv(
        args.out_file,
        index=False,
        columns=columns,
        quoting=csv.QUOTE_NONNUMERIC)
