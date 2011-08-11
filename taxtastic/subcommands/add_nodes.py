"""Add new nodes to a database containing a taxonomy."""

from taxtastic import ncbi
from taxtastic.taxonomy import Taxonomy
from taxtastic.utils import get_new_nodes

from sqlalchemy import create_engine
from sqlalchemy.exc import IntegrityError
import os.path

import logging
log = logging.getLogger(__name__)

def build_parser(parser):

    parser.add_argument(
        '-d', '--database-file',
        dest = 'database_file',
        metavar = 'FILE',
        required = True,
        help = 'Name of the sqlite database file')

    parser.add_argument(
        '-N', '--new-nodes',
        dest = 'new_nodes',
        metavar = 'FILE',
        required = True,
        help = """An Excel spreadsheet (.xls format only; requires
        xlrd) or csv-format file defining nodes to add to the
        taxonomy.  Mandatory fields include "tax_id", "parent_id",
        "rank", "tax_name"; optional fields include "source_name",
        "source_id". Other columns are ignored.""")

    parser.add_argument(
        '-S', '--source-name',
        dest='source_name',
        default='unknown',
        help="""Identifies the source for new nodes (will override
        "source_name" in the input provided by
        `--new-nodes`). [%(default)s]""")

def action(args):


    dbname = args.database_file
    new_nodes = args.new_nodes
    source_name = args.source_name

    engine = create_engine('sqlite:///%s' % dbname, echo=args.verbosity > 2)
    tax = Taxonomy(engine, ncbi.ranks)

    log.warning('adding new nodes')
    nodes = get_new_nodes(new_nodes)
    for d in nodes:
        if source_name:
            d['source_name'] = source_name
            try:
                tax.add_node(**d)
            except IntegrityError:
                log.info('node with tax_id %(tax_id)s already exists' % d)
            else:
                log.info('added new node with tax_id %(tax_id)s' % d)

    engine.dispose()
