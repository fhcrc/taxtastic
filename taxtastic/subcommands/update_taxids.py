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

Replaces tax_ids as specified in table 'merged' in the taxonomy
database. Use in preparation for ``taxit taxtable``. Takes sequence
info file as passed to ``taxit create --seq-info``

"""
import codecs
import csv
import logging
import sys

import sqlalchemy as sa

import taxtastic
from taxtastic.taxonomy import Taxonomy

log = logging.getLogger(__name__)


def build_parser(parser):
    parser.add_argument(
        'infile',
        type=taxtastic.utils.Opener('r'),
        help='Input file with taxids. Use "-" for stdin.')
    parser = taxtastic.utils.add_database_args(parser)
    parser.add_argument(
        '--delimiter',
        default=',',
        metavar='',
        type=lambda x: codecs.decode(str(x), 'unicode_escape'),
        help='Infile columns delimiter [%(default)s]')
    parser.add_argument(
        '--taxid-column',
        default='tax_id',
        metavar='',
        help='name of column or index if headerless containing '
             'tax_ids to be replaced [%(default)s]')
    parser.add_argument(
        '--unknowns',
        metavar='',
        type=taxtastic.utils.Opener('wt'),
        help='optional output file containing rows with unknown tax_ids '
             'having no replacements in merged table')
    parser.add_argument(
        '-a', '--unknown-action',
        choices=['drop', 'ignore', 'error'],
        default='error',
        help='action to perform for tax_ids with no replacement '
             'in merged table [%(default)s]')
    parser.add_argument(
        '-o', '--outfile',
        default=sys.stdout,
        metavar='',
        type=taxtastic.utils.Opener('wt'),
        help='Modified version of input file [stdout]')

    # not implemented for now
    # parser.add_argument(
    #     '--use-names', action='store_true', default=False,
    #     help='Use tax_name to assign replacement for unknown tax_ids'),
    # parser.add_argument(
    #     '--name-column', default='tax_name',
    #     help=('column to use for name lookup if --use-name '
    #           'is specified [%(default)s]'))


def action(args):
    header = next(args.infile).strip().split(args.delimiter)
    if args.taxid_column in header:
        taxid_column = header.index(args.taxid_column)
    elif len(header) == 1:
        taxid_column = 0
        header = None
        args.infile.seek(0)
    elif args.taxid_column.isnumeric():
        taxid_column = int(args.taxid_column) - 1
        header = None
        args.infile.seek(0)
    else:
        raise ValueError("No column " + args.taxid_column)

    drop = args.unknown_action == 'drop'
    error = args.unknown_action == 'error'
    ignore = args.unknown_action == 'ignore'

    # TODO: remove unless --use-names is implemented
    # if args.use_names and args.name_column not in header:
    #     raise ValueError("No column " + args.name_column)

    writer = csv.writer(args.outfile, delimiter=args.delimiter)
    if header:
        writer.writerow(header)

    if args.unknowns:
        unknowns = csv.writer(args.unknowns, delimiter=args.delimiter)
        if header:
            unknowns.writerow(header)

    engine = sa.create_engine(args.url, echo=args.verbosity > 3)
    tax = Taxonomy(engine, schema=args.schema)

    with tax.engine.connect() as con:
        log.info('reading table merged')
        q = 'select old_tax_id, new_tax_id from {merged}'.format(**tax.tables)
        result = con.execute(sa.text(q))
        mergedict = dict(result.fetchall())

        log.info('reading tax_ids from table nodes')
        result = con.execute(sa.text(
            'select tax_id from {nodes}'.format(**tax.tables)))
        all_tax_ids = {x[0] for x in result.fetchall()}

    log.info('reading input file')
    for row in csv.reader(args.infile, delimiter=args.delimiter):
        tax_id = row[taxid_column]

        if tax_id in all_tax_ids:
            pass  # write row without modification
        elif tax_id in mergedict:
            row[taxid_column] = mergedict[tax_id]
        else:  # tax_id is unknown
            if args.unknowns:
                unknowns.writerow(row)

            if ignore:
                pass
            elif drop:
                continue
            elif error:
                sys.exit(f'Error: tax_id {tax_id} is unknown')

        writer.writerow(row)
