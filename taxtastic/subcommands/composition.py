"""
Show taxonomic composition of a reference package.
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
from collections import Counter
import sys
import argparse

from Bio import SeqIO

from taxtastic import refpkg

log = logging.getLogger(__name__)

def build_parser(parser):
    parser.add_argument('refpkg', action='store', metavar='refpkg',
                        help='the reference package to operate on')
    parser.add_argument('-r', '--rank', default = 'species',
                        help = 'rank at which to show composition [%(default)s]')
    parser.add_argument('-o', '--outfile', default = sys.stdout, type = argparse.FileType('w'),
                        help = 'rank at which to show composition. Use --rank=tax_id to show original classifications [%(default)s]')


def action(args):
    log.info('loading reference package')

    pkg = refpkg.Refpkg(args.refpkg, create=False)

    with open(pkg.file_abspath('taxonomy')) as f:
        taxdict = {r['tax_id']:r for r in csv.DictReader(f)}

    counts = Counter()
    with open(pkg.file_abspath('seq_info')) as f:
        for row in csv.DictReader(f):
            # tax_id at the specified rank
            tax_id = taxdict[row['tax_id']][args.rank] if row['tax_id'] else ''
            tax_name = taxdict[tax_id]['tax_name'] if tax_id else '<unclassified at this rank>'
            counts[(tax_name, tax_id)] += 1

    writer = csv.writer(args.outfile)
    writer.writerow(['tax_name','tax_id','count'])
    writer.writerows(sorted((n, i, c) for (n, i),c in counts.items()))
