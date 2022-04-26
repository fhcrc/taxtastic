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
"""Find primary name and tax_id from taxonomic names

"""
import logging
import argparse
import sys
import csv

import sqlalchemy
import taxtastic

from taxtastic.taxonomy import Taxonomy

log = logging.getLogger(__name__)


def build_parser(parser):
    parser = taxtastic.utils.add_database_args(parser)

    input_group = parser.add_argument_group(
        "Input options").add_mutually_exclusive_group()

    input_group.add_argument(
        '-i', '--infile', metavar='FILE', type=argparse.FileType('r'),
        help='file containing a list of taxonomic names, one per line')

    input_group.add_argument(
        '-n', '--names', metavar='NAMES',
        help=('one or more taxonomic names provided as a comma-delimited '
              'list on the command line'))

    output_group = parser.add_argument_group("Output options")
    output_group.add_argument(
        '-o', '--outfile', metavar='FILE', type=argparse.FileType('w'),
        help='output file', default=sys.stdout)
    output_group.add_argument(
        '--include-unmatched', action='store_true', default=False,
        help='include names with no match')


def action(args):
    engine = sqlalchemy.create_engine(args.url, echo=False)
    tax = Taxonomy(engine, schema=args.schema)

    names = []
    if args.infile:
        names += [line.split('#', 1)[0].strip()
                  for line in args.infile
                  if line.strip() and not line.startswith('#')]

    if args.names:
        names += [x.strip() for x in args.names.split(',')]

    found_ids = dict(list(tax.id_from_names(names)))
    primary = dict(list(tax.primary_from_ids(found_ids.values())))

    writer = csv.writer(args.outfile)
    writer.writerow(['input', 'tax_name', 'tax_id', 'rank'])

    for name in names:
        if name in found_ids:
            tax_id = found_ids[name]
            _, rank = tax._node(tax_id)
            writer.writerow([name, primary[tax_id], tax_id, rank])
        elif args.include_unmatched:
            writer.writerow([name, None, None, None])
        else:
            log.warning('dropping ({}), not found in database'.format(name))

    log.warning('found {} of {} names'.format(len(found_ids), len(names)))
