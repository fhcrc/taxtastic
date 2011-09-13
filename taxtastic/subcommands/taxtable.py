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

import logging
import argparse

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

    dbname = args.database_file
    taxids = args.taxids
    taxnames = args.taxnames
    csvfile = args.out_file

    engine = create_engine('sqlite:///%s' % dbname, echo=args.verbosity > 2)
    tax = Taxonomy(engine, ncbi.ranks)

    # get a list of taxa
    taxa = set()

    if taxids:
        if os.access(taxids, os.F_OK):
            log.warning('reading tax_ids from %s' % taxids)
            for line in getlines(taxids):
                taxa.update(set(line.split()))
        else:
            taxa = set([x.strip() for x in taxids.split(',')])

    if taxnames:
        for tax_name in getlines(taxnames):
            tax_id, primary_name, is_primary = tax.primary_from_name(tax_name)
            taxa.add(tax_id)
            if not is_primary:
                log.warning(
                    '%(tax_id)8s  %(tax_name)40s -(primary name)-> %(primary_name)s' \
                        % locals())

    log.warning('calculating lineages for %s taxa' % len(taxa))
    for taxid in taxa:
        log.info('adding %s' % taxid)
        tax.lineage(taxid)

    log.warning('writing output to %s' % csvfile.name)
    tax.write_table(None, csvfile = csvfile)

    engine.dispose()
