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
"""Add new nodes to a database"""
from taxtastic.taxonomy import Taxonomy, TaxonIntegrityError
from taxtastic.utils import get_new_nodes, add_database_args

import logging
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

    parser = add_database_args(parser)

    parser.add_argument(
        '-S', '--source-name',
        dest='source_name',
        default='unknown',
        help=('Identifies the source for new nodes (will override '
              '"source_name" in the input provided '
              'by `--new-nodes`). [%(default)s]'))

    parser.add_argument(
        '--update',
        action='store_true',
        help='Update any existing nodes')

    parser.add_argument(
        '--tax-table',
        help='table name in database')

    parser.add_argument(
        '--out',
        default=sys.stdout,
        help='sql queries')


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
    engine = sqlalchemy.create_engine(args.url, echo=args.verbosity > 2)
    tax = Taxonomy(engine, schema=args.schema)
    nodes = list(get_new_nodes(args.new_nodes))

    # check if there are any new ranks and exit if needed
    node_ranks = set(n['rank'] for n in nodes)
    for r in node_ranks:
        if r not in tax.ranks:
            msg = 'adding new ranks to taxonomy is not yet supported'
            raise TaxonIntegrityError(msg)

    ranksdict = tax.ranksdict()
    ranksdict.update(dict([(n['tax_id'], n['rank']) for n in nodes]))
    nodes = [verify_rank_integrity(n, ranksdict, tax.ranks) for n in nodes]
    nodes = [verify_lineage_integrity(n, ranksdict, tax.ranks, tax) for n in nodes]

    log.info('adding new nodes')
    for d in nodes:
        if args.source_name:
            d['source_name'] = args.source_name
            try:
                tax.add_node(**d)
            except sqlalchemy.exc.IntegrityError:
                if args.update:
                    tax.update_node(**d)
                else:
                    log.warn('node with tax_id %(tax_id)s already exists' % d)
            else:
                log.info('added new node with tax_id %(tax_id)s' % d)

    engine.dispose()
