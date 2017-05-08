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
import taxtastic

log = logging.getLogger(__name__)


def build_parser(parser):
    parser.add_argument(
        'infile',
        nargs='?',
        default=sys.stdin,
        help=('Input CSV file to process, minimally '
              'containing the fields `tax_id`. Rows with '
              'missing tax_ids are left unchanged.'))

    parser = taxtastic.utils.add_database_args(parser)

    parser.add_argument(
        '-o', '--out',
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


def action(args):
    try:
        rows = pandas.read_csv(args.infile, dtype='str')
    except pandas.io.common.EmptyDataError as e:
        log.error(e)
        return
    columns = rows.columns.tolist()  # preserve column order

    if args.taxid_column not in columns:
        raise ValueError("No column " + args.taxid_column)

    if args.name_column:
        if args.name_column not in columns:
            msg = '"No "' + args.name_column + '" column'
            raise ValueError(msg)

    engine = sqlalchemy.create_engine(args.url, echo=args.verbosity > 3)

    merged = pandas.read_sql_table(
        'merged', engine,
        schema=args.schema,
        index_col='old_tax_id')
    log.info('updating tax_ids')
    rows = rows.join(merged, on=args.taxid_column)

    # overwrite tax_ids where there is a new_tax_id
    inew_tax_ids = ~rows['new_tax_id'].isnull()
    rows.loc[inew_tax_ids, args.taxid_column] = \
        rows[inew_tax_ids]['new_tax_id']
    rows = rows.drop('new_tax_id', axis=1)

    log.info('loading names table')
    names = pandas.read_sql_table(
        'names', engine,
        schema=args.schema,
        columns=['tax_id', 'tax_name', 'is_primary'])

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
        split seq_info into two dfs rank <> species

        ranks <= species get is_valid column
        ranks > species get False
        """

        if 'taxid_classified' in columns:
            rows = rows.drop('taxid_classified', axis=1)
        else:
            columns.append('taxid_classified')
        ranks = pandas.read_sql_table('ranks', engine, schema=args.schema)
        nodes = pandas.read_sql_table(
            'nodes', engine,
            schema=args.schema,
            columns=['tax_id', 'rank', 'is_valid'],
            index_col='tax_id')
        rows = species_is_classified(rows, nodes, ranks)

    # output seq_info with new tax_ids
    rows.to_csv(
        args.out,
        index=False,
        columns=columns,
        quoting=csv.QUOTE_NONNUMERIC)


def species_is_classified(rows, nodes, ranks):
    rows = rows.join(nodes, on='tax_id')
    rows['rank'] = rows['rank'].astype(
        'category', categories=ranks['rank'].tolist(), ordered=True)
    rows.loc[rows['rank'] > 'species', 'is_valid'] = False
    rows = rows.rename(columns={'is_valid': 'taxid_classified'})
    return rows
