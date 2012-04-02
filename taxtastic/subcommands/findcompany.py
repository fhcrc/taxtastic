"""Find company for lonely nodes."""
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
import logging

from taxtastic import lonely
from sqlalchemy import create_engine
from taxtastic.taxonomy import Taxonomy
from taxtastic import ncbi


log = logging.getLogger(__name__)

class ConfigError(Exception):
    pass

def build_parser(parser):
    parser.add_argument("taxdb",
                        help="Taxonomy database to work from")
    parser.add_argument("tax_ids", type=str, nargs='*',
                        help='Tax IDs to look up')
    parser.add_argument("-c", "--cut",
                        action="store_true", default=False,
                        help="Return one output per input (probably also lonely)")
    parser.add_argument("-i", "--input", type=argparse.FileType('r'),
            default=None, help="Text file to read Tax IDs from, one per line")
    parser.add_argument('-o', '--output',
                        default=None,
                        help='Write output to given file')


def action(args):
    taxids = args.tax_ids
    # Add taxids from input file
    if args.input:
        with args.input as h:
            for l in h:
                val = l.split('#')[0].strip()
                taxids.append(val)
    # Connect to the taxonomy
    engine = create_engine('sqlite:///%s' % args.taxdb, echo=False)
    tax = Taxonomy(engine, ncbi.ranks)
    # Finally, real work...
    if args.cut:
        company = lonely.lonely_company(tax, taxids)
    else:
        company = lonely.solid_company(tax, taxids)
    txt = ""
    for t in company:
        txt += "%s\n" % (t if t else "")
    if args.output:
        with open(args.output, 'w') as h:
            print >>h, txt
    else:
        print txt
    return 0
