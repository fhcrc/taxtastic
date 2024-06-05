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
"""
Returns given taxids including descendant taxids
"""
import argparse
import logging
import os
import sqlalchemy
import sys
import taxtastic

from taxtastic.taxonomy import Taxonomy

log = logging.getLogger(__name__)


def build_parser(parser):
    parser = taxtastic.utils.add_database_args(parser)
    parser.add_argument(
        'taxids',
        help='File or comma delimited list of tax_ids')
    parser.add_argument(
        '--out',
        type=argparse.FileType('w'),
        default=sys.stdout)


def action(args):
    if os.path.isfile(args.taxids) or args.taxids == '-':
        taxids = [i.strip() for i in taxtastic.utils.Opener('rt')(args.taxids) if i]
    else:
        taxids = args.taxids.split(',')
    engine = sqlalchemy.create_engine(args.url, echo=args.verbosity > 3)
    tax = Taxonomy(engine, schema=args.schema)
    for i in tax.descendants_of(taxids):
        args.out.write(i + '\n')
