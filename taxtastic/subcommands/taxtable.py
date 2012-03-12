"""Creates a CSV file describing lineages for a set of taxa"""
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

import re

from taxtastic import ncbi
from taxtastic.taxonomy import Taxonomy
from taxtastic.utils import getlines

from sqlalchemy import create_engine
from sqlalchemy.exc import IntegrityError
import os.path
import sys

log = logging.getLogger(__name__)

def build_parser(parser):

    parser.add_argument(
        '-d', '--database-file',
        dest = 'database_file',
        metavar = 'FILE',
        required = True,
        help = 'Name of the sqlite database file')

    input_group = parser.add_argument_group(
        "Input options").add_mutually_exclusive_group()

    input_group.add_argument(
        '-n', '--tax-names',
        dest = 'taxnames',
        metavar = 'FILE',
        help="""A file identifing taxa in the form of taxonomic
        names. Names are matched against both primary names and
        synonyms. Lines beginning with "#" are ignored. Taxa
        identified here will be added to those specified using
        --tax-ids""")

    input_group.add_argument(
        '-t', '--tax-ids',
        dest = 'taxids',
        metavar = 'FILE-OR-LIST',
        help = """File containing a whitespace-delimited list of
        tax_ids (ie, separated by tabs, spaces, or newlines; lines
        beginning with "#" are ignored). This option can also be
        passed a comma-delited list of taxids on the command line.""")

    input_group.add_argument(
        '-i', '--seq-info', type=argparse.FileType('r'),
        help="""Read tax_ids from sequence info file, minimally containing the
        field "tax_id" """)

    output_group = parser.add_argument_group(
        "Output options").add_mutually_exclusive_group()

    output_group.add_argument(
        '-o', '--out-file',
        dest='out_file',
        type = argparse.FileType('w'),
        default = sys.stdout,
        metavar='FILE',
        help="""Output file containing lineages for the specified taxa
        in csv format; writes to stdout if unspecified""")

def action(args):
    engine = create_engine('sqlite:///%s' % args.database_file, echo=args.verbosity > 2)
    tax = Taxonomy(engine, ncbi.ranks)

    taxids = set()

    if args.taxids:
        if os.access(args.taxids, os.F_OK):
            for line in getlines(args.taxids):
                taxids.update(set(re.split(r'[\s,;]+', line)))
        else:
            taxids.update([x.strip() for x in re.split(r'[\s,;]+', args.taxids)])

    if args.taxnames:
        for taxname in getlines(args.taxnames):
            for name in re.split(r'\s*[,;]\s*', taxname):
                tax_id, primary_name, is_primary = tax.primary_from_name(name.strip())
                taxids.add(tax_id)

    if args.seq_info:
        with args.seq_info:
            reader = csv.DictReader(args.seq_info)
            taxids.update(frozenset(i['tax_id'] for i in reader if i['tax_id']))

    # Before digging into lineages, make sure all the taxids exist in
    # the taxonomy database.
    valid_taxids = True
    for t in taxids:
        try:
            tax._node(t)
        except KeyError:
            # Check for merged
            m = tax._get_merged(t)
            if m and m != t:
                print >> sys.stderr, ("Taxid {0} has been replaced by {1}. "
                        "Please update your records").format(t, m)
            else:
                print >>sys.stderr, "Taxid %s not found in taxonomy." % t
            valid_taxids = False
    if not(valid_taxids):
        print >>sys.stderr, "Some taxids were invalid.  Exiting."
        return 1 # exits with code 1

    # Extract all the taxids to be exported in the CSV file.
    taxids_to_export = set()
    for t in taxids:
        taxids_to_export.update([y for (x,y) in tax._get_lineage(t)])

    tax.write_table(taxids_to_export, csvfile = args.out_file)

    engine.dispose()
    return 0
