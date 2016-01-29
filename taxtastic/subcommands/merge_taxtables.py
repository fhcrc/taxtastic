"""Combine two or more taxtables

The first file listed identifies the ranks to be included in the
output. In subsequent files, the introduction of additional ranks not
found in the first file causes an error.

"""

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

import argparse
import csv
import logging
import sys
import itertools

log = logging.getLogger(__name__)


def build_parser(parser):
    parser.add_argument(
        'infiles', nargs='+', type=argparse.FileType(),
        help="""One or more taxtables""")
    parser.add_argument(
        '-o', '--outfile', type=argparse.FileType('w'), default=sys.stdout,
        metavar='FILE', help="""csv file containing merged taxtables""")


def action(args):
    first, others = args.infiles[0], args.infiles[1:]
    if not others:
        sys.exit('At least two taxtables are required')

    reader = csv.DictReader(first)
    fieldnames = reader.fieldnames[:]
    initial_names = set(fieldnames)
    writer = csv.DictWriter(args.outfile, fieldnames=reader.fieldnames)
    writer.writeheader()
    empty_row = {f: None for f in reader.fieldnames}

    tax_ids = {}
    for lineage in itertools.chain(reader, *(csv.DictReader(f) for f in others)):
        if set(lineage.keys()) - initial_names:
            log.error('A file contains taxa not found in {}: {}'.format(
                first.name, set(lineage.keys()) - initial_names))
            sys.exit(1)

        tax_id = lineage['tax_id']
        new_lineage = dict(empty_row, **{k: v or None for k, v in lineage.items()})
        if tax_id in tax_ids:
            existing_lineage = tax_ids[tax_id]
            if existing_lineage != new_lineage:
                log.error('lineages for {} do not match'.format(tax_id))
                for key in fieldnames:
                    v1, v2 = existing_lineage[key], new_lineage[key]
                    if v1 != v2:
                        print [key, v1, v2]
        else:
            writer.writerow(new_lineage)
            tax_ids[tax_id] = new_lineage

