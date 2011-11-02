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
import os.path
import shutil
import csv
from collections import defaultdict

from Bio import Phylo, SeqIO

from taxtastic import refpkg

log = logging.getLogger(__name__)

def build_parser(parser):
    parser.add_argument('refpkg', action='store', metavar='refpkg',
                        help='the reference package to operate on')
    parser.add_argument('-n', '--seq-names', action = 'store_true', default = False,
                        help = 'print a list of sequence names')
    parser.add_argument('-t', '--taxonomy', action = 'store_true', default = False,
                        help = 'print a tally of sequences representing each taxon at rank RANK')
    parser.add_argument('-r', '--rank', metavar = 'RANK', default = 'species',
                        help = 'rank at which to describe the taxonomy [%(default)s]')

def describe_taxonomy(pkg, rank):
    tally = defaultdict(int)
    with open(pkg.file_abspath('taxonomy')) as taxtab, open(pkg.file_abspath('seq_info')) as seq_info:  
        taxdict = {row['tax_id']: row for row in csv.DictReader(taxtab)}
        taxdict[''] = {'tax_name':'undefined'}
        for refseq in csv.DictReader(seq_info):
            tally[taxdict[refseq['tax_id']][rank]] += 1

    items = [(taxdict[tax_id]['tax_name'], count) for tax_id, count in tally.items()]
    for name, count in sorted(items):
        print name, count
        
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
    elif args.taxonomy:
        describe_taxonomy(pkg, args.rank)
    else:
        print 'number of sequences:', len(snames)
        print 'package components\n', '\n'.join(sorted(pkg.file_keys())) 
