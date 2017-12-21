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
"""Create a table of lineages as taxonimic names for a collection of sequences

Minimal inputs are a taxtable and a file providing a mapping of
sequence names to tax_ids. Outputs are one or more of:

* a table of taxonomic lineges in csv format

* a "taxonomy" file formatted for MOTHUR
  (https://mothur.org/wiki/Taxonomy_File). Ranks are limited to the
  following, with corresponding abbreviations:

  [('superkingdom', 'pk'),
   ('phylum', 'ph'),
   ('class', 'cl'),
   ('order', 'or'),
   ('family', 'fa'),
   ('genus', 'ge'),
   ('species', 'sp')]

  Lineages are truncated to species, and missing values reported
  as '<abbrev>__unclassified'

"""

import argparse
import csv
import logging
from operator import itemgetter
from collections import OrderedDict

log = logging.getLogger(__name__)


def build_parser(parser):
    input_group = parser.add_argument_group('input options')

    input_group.add_argument(
        'taxtable', metavar='FILE', type=argparse.FileType('rt'),
        help=('output of "taxit taxtable" containing all tax_ids represented in "seq_info"'))

    input_group.add_argument(
        'seq_info', metavar='FILE',
        type=argparse.FileType('rt'),
        help=('csv file providing a mapping of sequence names to tax_ids '))

    input_group.add_argument(
        '--seqname-col', metavar='NAME', default='seqname',
        help=('name of column in "seq_info" containing sequence names'))

    input_group.add_argument(
        '--tax-id-col', metavar='NAME', default='tax_id',
        help=('name of column in "seq_info" containing tax_ids'))

    output_group = parser.add_argument_group(
        "Output options").add_mutually_exclusive_group()

    output_group.add_argument(
        '-c', '--csv-table', type=argparse.FileType('wt'), metavar='FILE',
        help=('Output file containing lineages for each sequence name '
              'in csv format'))

    output_group.add_argument(
        '-t', '--taxonomy-table', type=argparse.FileType('wt'), metavar='FILE',
        help=('"taxonomy" file formatted for MOTHUR'))


def action(args):

    # mapping of seqname: tax_id in input order
    namedict = OrderedDict(map(
        itemgetter(args.seqname_col, args.tax_id_col), csv.DictReader(args.seq_info)))

    log.info('read {} sequence names'.format(len(namedict)))

    reader = csv.DictReader(args.taxtable)
    taxdict = {row['tax_id']: row for row in reader}
    taxnames = {d['tax_id']: d['tax_name'] for d in taxdict.values()}

    taxcols = reader.fieldnames
    all_ranks = taxcols[taxcols.index('root'):]

    lineages = {}
    for tax_id in set(namedict.values()):
        d = taxdict[tax_id]
        lineages[tax_id] = {rank: taxnames.get(d[rank]) for rank in all_ranks}

    if args.csv_table:
        writer = csv.DictWriter(args.csv_table, fieldnames=['seqname'] + all_ranks)
        writer.writeheader()
        for name, tax_id in namedict.items():
            writer.writerow(dict(lineages[tax_id], seqname=name))

    if args.taxonomy_table:
        # abbreviations from http://blog.mothur.org/2017/03/22/SILVA-v128-reference-files/
        base_ranks = [('superkingdom', 'pk'),
                      ('phylum', 'ph'),
                      ('class', 'cl'),
                      ('order', 'or'),
                      ('family', 'fa'),
                      ('genus', 'ge'),
                      ('species', 'sp')]

        for name, tax_id in namedict.items():
            lineage = lineages[tax_id]
            taxstr = ';'.join('"{}__{}"'.format(abbrev, lineage.get(rank) or 'unclassified')
                              for rank, abbrev in base_ranks)
            args.taxonomy_table.write('{}\t{}\n'.format(name, taxstr))
