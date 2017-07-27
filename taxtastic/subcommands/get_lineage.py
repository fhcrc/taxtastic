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
"""Calculate the taxonomic lineage of a taxid
"""

import argparse
import csv
import logging
import sys
from itertools import islice

import sqlalchemy

from taxtastic.taxonomy import Taxonomy
from taxtastic.utils import add_database_args

log = logging.getLogger(__name__)


def build_parser(parser):
    parser = add_database_args(parser)
    # parser.add_argument('tax_id')
    parser.add_argument('tax_ids', type=argparse.FileType())
    parser.add_argument(
        '-o', '--outfile',
        type=argparse.FileType('w'),
        default=sys.stdout,
        metavar='FILE')


def action(args):

    engine = sqlalchemy.create_engine(args.url, echo=args.verbosity > 3)
    tax = Taxonomy(engine, schema=args.schema)

    for tax_id in islice(args.tax_ids, 5):
        print tax._get_lineage_table(tax_id.strip())

    # print tax._get_lineage(args.tax_id)
    # print tax.lineage(args.tax_id)

    engine.dispose()
