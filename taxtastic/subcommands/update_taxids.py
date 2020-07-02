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

import csv
import logging
import sqlalchemy
import sys
import taxtastic

from fastalite import Opener

from taxtastic.taxonomy import Taxonomy

log = logging.getLogger(__name__)


def build_parser(parser):
    parser.add_argument(
        'infile', type=Opener('rU'),
        help=('Input CSV file to process, minimally containing the field `tax_id`. '
              'Use "-" for stdin.'))
    parser = taxtastic.utils.add_database_args(parser)
    parser.add_argument(
        '-o', '--outfile', default=sys.stdout, type=Opener('wt'),
        help='Modified version of input file [default: stdout]')
    parser.add_argument(
        '--taxid-column', default='tax_id',
        help='name of column containing tax_ids to be replaced [%(default)s]')
    parser.add_argument(
        '--unknowns', type=Opener('wt'),
        help=('optional output file containing rows with unknown tax_ids '
              'having no replacements in merged table'))
    parser.add_argument(
        '-a', '--unknown-action', choices=['drop', 'ignore', 'error'], default='error',
        help=('action to perform for tax_ids with no replacement '
              'in merged table [%(default)s]'))

    # not implemented for now
    # parser.add_argument(
    #     '--use-names', action='store_true', default=False,
    #     help='Use tax_name to assign replacement for unknown tax_ids'),
    # parser.add_argument(
    #     '--name-column', default='tax_name',
    #     help=('column to use for name lookup if --use-name '
    #           'is specified [%(default)s]'))


def action(args):
    reader = csv.DictReader(args.infile)
    fieldnames = reader.fieldnames
    taxid_column = args.taxid_column
    drop = args.unknown_action == 'drop'
    error = args.unknown_action == 'error'
    ignore = args.unknown_action == 'ignore'

    if taxid_column not in fieldnames:
        raise ValueError("No column " + args.taxid_column)

    # TODO: remove unless --use-names is implemented
    # if args.use_names and args.name_column not in fieldnames:
    #     raise ValueError("No column " + args.name_column)

    writer = csv.DictWriter(args.outfile, fieldnames=fieldnames, quoting=csv.QUOTE_ALL)
    writer.writeheader()

    if args.unknowns:
        unknowns = csv.DictWriter(
            args.unknowns, fieldnames=fieldnames, quoting=csv.QUOTE_ALL)
        unknowns.writeheader()

    engine = sqlalchemy.create_engine(args.url, echo=args.verbosity > 3)
    tax = Taxonomy(engine, schema=args.schema)

    with tax.engine.connect() as con:
        log.info('reading table merged')
        result = con.execute(
            'select old_tax_id, new_tax_id from {merged}'.format(merged=tax.merged))
        mergedict = dict(result.fetchall())

        log.info('reading tax_ids from table {nodes}'.format(nodes=tax.nodes))
        result = con.execute('select tax_id from {nodes}'.format(nodes=tax.nodes))
        all_tax_ids = {x[0] for x in result.fetchall()}

    log.info('reading input file')
    for row in reader:
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
                sys.exit('Error: tax_id {} is unknown'.format(tax_id))

        writer.writerow(row)
