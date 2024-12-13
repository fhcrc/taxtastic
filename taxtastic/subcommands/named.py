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
Filters unclassified, unnamed taxonomy ids
"""
import csv
import sqlalchemy
import sys
from taxtastic.utils import add_database_args, Opener
from taxtastic.taxonomy import Taxonomy


def build_parser(parser):
    parser = add_database_args(parser)
    input_group = parser.add_mutually_exclusive_group(required=True)
    input_group.add_argument(
        '-t', '--tax-ids',
        nargs='+',
        help='one or more space-delimited tax_ids (eg "-t 47770 33945")')
    input_group.add_argument(
        '-f', '--tax-id-file',
        metavar='FILE',
        type=Opener('rt'),
        help=('File containing a whitespace-delimited list of '
              'tax_ids (ie, separated by tabs, spaces, or newlines.'))
    input_group.add_argument(
        '-i', '--seq-info',
        type=Opener('rt'),
        help=('Read tax_ids from sequence info file, minimally '
              'containing a column named "tax_id"'))
    parser.add_argument(
        '--ranked',
        action='store_true',
        help='Ignore "no rank" taxonomies [%(default)s]')
    parser.add_argument(
        '-o', '--outfile',
        type=Opener('wt'),
        default=sys.stdout,
        metavar='FILE',
        help=('Output file containing named taxonomy ids;'
              'writes to stdout if unspecified'))


def action(args):
    engine = sqlalchemy.create_engine(args.url, echo=args.verbosity > 3)
    tax = Taxonomy(engine, schema=args.schema)
    if args.seq_info:
        seq_info = csv.DictReader(args.seq_info)
        fieldnames = seq_info.fieldnames
        seq_info = list(seq_info)
        tax_ids = (i['tax_id'] for i in seq_info)
        named = set(tax.is_valid(tax_ids, no_rank=not args.ranked))
        out = csv.DictWriter(args.outfile, fieldnames=fieldnames)
        out.writeheader()
        out.writerows(i for i in seq_info if i['tax_id'] in named)
    else:
        if args.tax_ids:
            tax_ids = args.tax_ids
        else:
            tax_ids = (i.strip() for i in args.tax_id_file)
            tax_ids = [i for i in tax_ids if i]
        named = set(tax.is_valid(tax_ids, no_rank=not args.ranked))
        tax_ids = (i for i in tax_ids if i in named)
        for i in tax_ids:
            args.outfile.write(i + '\n')
