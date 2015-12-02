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
        log.info('name {} associated with tax_id {}'.format(name, tax_id))
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

    tax_ids = rows[['tax_id']].drop_duplicates()
    halt = not (bool(args.unknowns) or bool(args.name_column))

    def set_new_tax_id(row):
        row['new_tax_id'] = update_taxid(row['tax_id'], tax, halt=halt)
        return row

    msg = 'updating taxids:'
    tax_ids = utils.apply_df_status(set_new_tax_id, tax_ids, msg)

    if args.name_column:
        """
        Seperate new_tax_ids that are not found and use the name_column
        to search for an id.
        """
        null_ids = tax_ids['new_tax_id'].isnull()
        no_tax_ids = tax_ids[null_ids]
        no_tax_ids = no_tax_ids.drop('new_tax_id', axis=1)  # drop empty column
        all_names = rows[['tax_id', args.name_column]].drop_duplicates()
        no_tax_ids = no_tax_ids.merge(all_names, on='tax_id', how='left')
        names = no_tax_ids[[args.name_column]].drop_duplicates()

        halt = not bool(args.unknowns)

        def set_name_taxid(row):
            row['new_tax_id'] = update_by_name(
                row[args.name_column], tax, halt=halt)
            return row

        msg = 'updating tax_ids by ' + args.name_column
        names = utils.apply_df_status(set_name_taxid, names, msg)

        no_tax_ids = no_tax_ids.merge(names, on=args.name_column, how='left')
        no_tax_ids = no_tax_ids.drop(args.name_column, axis=1)
        tax_ids = [tax_ids[~null_ids], no_tax_ids]
        tax_ids = pandas.concat(tax_ids)

    # remove new_tax_id.isnull()
    tax_ids = tax_ids[~tax_ids['new_tax_id'].isnull()]

    new_tax_ids = tax_ids[['new_tax_id']].drop_duplicates()
    if args.taxid_classified:
        """
        Decide if tax_id is a valid taxon
        """
        if 'taxid_classified' in columns:
            rows = rows.drop('taxid_classified', axis=1)
        else:
            columns.append('taxid_classified')

        new_tax_ids.loc[:, 'taxid_classified'] = None

        def is_classified(row):
            row['taxid_classified'] = species_is_classified(
                row['new_tax_id'], tax)
            return row

        msg = 'validating tax_ids:'
        new_tax_ids = utils.apply_df_status(is_classified, new_tax_ids, msg)

    if args.append_lineage:
        """
        Append a column from the taxonomy to seq_info
        """
        if args.append_lineage in columns:
            rows = rows.drop(args.append_lineage, axis=1)
        else:
            columns.append(args.append_lineage)

        new_tax_ids.loc[:, args.append_lineage] = None

        def add_rank_column(row):
            lineage = tax.lineage(row['new_tax_id'])
            row[args.append_lineage] = lineage.get(args.append_lineage, None)
            return row

        msg = 'appending {} column'.format(args.append_lineage)
        new_tax_ids = utils.apply_df_status(add_rank_column, new_tax_ids, msg)

    tax_ids = tax_ids.merge(new_tax_ids, on='new_tax_id', how='left')
    rows = rows.merge(tax_ids, on='tax_id', how='left')

    unknowns = rows['new_tax_id'].isnull()
    if args.unknowns:
        unknown_rows = rows[unknowns].drop('new_tax_id', axis=1)
        unknown_rows.to_csv(args.unknowns,
                            index=False, columns=columns,
                            quoting=csv.QUOTE_NONNUMERIC)

    # output seq_info with new tax_ids
    known_rows = rows[~unknowns].drop('tax_id', axis=1)
    known_rows = known_rows.rename(columns={'new_tax_id': 'tax_id'})
    known_rows.to_csv(args.out_file,
                      index=False, columns=columns,
                      quoting=csv.QUOTE_NONNUMERIC)
