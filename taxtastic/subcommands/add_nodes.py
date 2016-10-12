"""Add new nodes to a database containing a taxonomy."""

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

from taxtastic.taxonomy import Taxonomy, TaxonIntegrityError
from taxtastic.utils import get_new_nodes

import logging
import pandas
import sqlalchemy
import sys

log = logging.getLogger(__name__)


def build_parser(parser):

    parser.add_argument(
        'new_nodes',
        metavar='csv',
        help=('A csv file defining nodes to add to the taxonomy. '
              'Mandatory fields include "tax_id", "parent_id", "rank", '
              '"tax_name"; optional fields include "source_name", '
              '"source_id" and "children". The "children" field should '
              'specify one or more existing taxids in a semicolon-'
              'delimited list. Other columns are ignored.'))

    parser.add_argument(
        'database_file',
        metavar='sqlite',
        help='Name of the database file')

    parser.add_argument(
        '--no-write-to-database-file',
        dest='write_to_database_file',
        action='store_false',
        help=('do not results to database file'))

    parser.add_argument(
        '-S', '--source-name',
        dest='source_name',
        default='unknown',
        help=('Identifies the source for new nodes (will override '
              '"source_name" in the input provided '
              'by `--new-nodes`). [%(default)s]'))

    table_parser = parser.add_argument_group('taxtable')
    table_parser.add_argument(
        '--from-table',
        metavar='csv',
        help=('a taxonomic table'))
    table_parser.add_argument(
        '--to-table',
        metavar='csv',
        default=sys.stdout,
        help=('output taxonomic table.  requires --from-table'))


def verify_rank_integrity(node, ranksdict, rank_order):
    '''
    confirm that for each node the parent ranks and children ranks are coherent
    '''
    def _lower(n1, n2):
        return rank_order.index(n1) < rank_order.index(n2)

    rank = node['rank']
    if not _lower(rank, ranksdict[node['parent_id']]):
        msg = 'New node {} has same or higher rank than parent node {}'
        msg = msg.format(node['tax_id'], node['parent_id'])
        raise TaxonIntegrityError(msg)
    if 'children' in node:
        for child in node['children']:
            if not _lower(ranksdict[child], rank):
                msg = 'Child node {} has same or lower rank as new node {}'
                msg = msg.format(node['tax_id'], child)
                raise TaxonIntegrityError(msg)

    return node


def verify_lineage_integrity(node, ranksdict, rank_order, tax):
    rank = node['rank']

    if 'children' in node:
        def _le(r1, r2):
            return rank_order.index(r1) <= rank_order.index(r2)

        for c in node['children']:
            try:
                parent_parent_id = tax.parent_id(tax.parent_id(c))
            except ValueError:
                continue
            parent_parent_rank = ranksdict[parent_parent_id]
            if _le(parent_parent_rank, rank):
                msg = 'New node {} splits lineage of existing child node {}'
                msg = msg.format(node['tax_id'], c)
                raise TaxonIntegrityError(msg)

    return node


def action(args):
    engine = sqlalchemy.create_engine(
        'sqlite:///' + args.database_file, echo=args.verbosity > 2)
    ranks = pandas.read_sql_table('ranks', engine)['rank'].tolist()
    nodes = list(get_new_nodes(args.new_nodes))

    # check if there are any new ranks and exit if needed
    node_ranks = set(n['rank'] for n in nodes)
    for r in node_ranks:
        if r not in ranks:
            msg = 'adding new ranks to taxonomy is not yet supported'
            raise TaxonIntegrityError(msg)

    tax = Taxonomy(engine, ranks)
    ranksdict = tax.ranksdict()
    nodes = [verify_rank_integrity(n, ranksdict, ranks) for n in nodes]
    nodes = [verify_lineage_integrity(n, ranksdict, ranks, tax) for n in nodes]

    if args.write_to_database_file:
        log.warn('adding new nodes')
        for d in nodes:
            if args.source_name:
                d['source_name'] = args.source_name
                try:
                    tax.add_node(**d)
                except sqlalchemy.exc.IntegrityError:
                    log.warn('node with tax_id %(tax_id)s already exists' % d)
                else:
                    log.info('added new node with tax_id %(tax_id)s' % d)

    if args.from_table:
        log.info('reading from_table ' + args.from_table)
        table = pandas.read_csv(args.from_table, dtype=str).set_index('tax_id')
        table['rank'] = table['rank'].astype('category', categories=ranks)
        from_cols = table.columns  # preserve order of columns
        nodes_df = pandas.DataFrame(nodes, dtype=str).set_index('tax_id')
        lineage_cols = [c for c in table.columns if c in ranks]

        log.info('appending new nodes to from_table ' + args.from_table)
        while not nodes_df.empty:
            new_rows = nodes_df.join(
                table[lineage_cols], how='inner', on='parent_id')
            for i, r in new_rows.iterrows():
                new_rows.loc[i, r['rank']] = i
            table = table.append(new_rows)
            nodes_df = nodes_df.drop(new_rows.index)

        log.info('updating new node children (of children)')
        for n in nodes:
            if 'children' in n:
                rank = n['rank']
                tax_id = n['tax_id']
                for c in n['children']:
                    c_rank = table.loc[c]['rank']
                    table.loc[table[c_rank] == c, rank] = tax_id

        table = table.sort_values(by='rank', ascending=False)
        table.to_csv(args.to_table, columns=from_cols)

    engine.dispose()
