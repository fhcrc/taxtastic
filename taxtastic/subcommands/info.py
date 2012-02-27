"""
Show information about reference packages.
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

import logging
import csv
from collections import defaultdict
import pprint
import sys

from Bio import Phylo

from taxtastic import refpkg

log = logging.getLogger(__name__)

def build_parser(parser):
    parser.add_argument('refpkg', action='store', metavar='refpkg',
                        help='the reference package to operate on')
    parser.add_argument('-n', '--seq-names', action = 'store_true', default = False,
                        help = 'print a list of sequence names')
    parser.add_argument('-t', '--tally', action = 'store_true', default = False,
                        help = 'print a tally of sequences representing each taxon at rank RANK')

def tally_taxa(pkg):
    tally = defaultdict(int)
    with open(pkg.file_abspath('taxonomy')) as taxtab, open(pkg.file_abspath('seq_info')) as seq_info:
        taxdict = {row['tax_id']: row for row in csv.DictReader(taxtab)}

        tax_ids = [d['tax_id'] for d in csv.DictReader(seq_info)]

    for tax_id in tax_ids:
        tally[tax_id] += 1
        
    rows = [(taxdict[tax_id]['tax_name'], tax_id, count) for tax_id, count in tally.items()]

    writer = csv.writer(sys.stdout, quoting = csv.QUOTE_NONNUMERIC)
    writer.writerows(sorted(rows))

def action(args):
    """
    Show information about reference packages.
    """
    log.info('loading reference package')

    pkg = refpkg.Refpkg(args.refpkg)

    with open(pkg.file_abspath('seq_info')) as seq_info:
        snames = [row['seqname'] for row in csv.DictReader(seq_info)]

    if args.seq_names:
        print '\n'.join(snames)
    elif args.tally:
        tally_taxa(pkg)
    else:
        print 'number of sequences:', len(snames)
        print 'package components\n', '\n'.join(sorted(pkg.file_keys()))
