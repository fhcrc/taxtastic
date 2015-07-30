"""Look up tax_ids from ncbi from accession numbers"""
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

import csv
import logging
import argparse
import sys

log = logging.getLogger(__name__)


def id_from_accessions(accessions, email=None):
    """
    Retrieve tax ids from a list of accessions.

    ncbi.efetch can work unexpectedly when querying a large number of tax ids.
    The best way to mitigate is to send only its documented max of 10k ids
    at a time.  Ncbi claims it can handle more than 10k ids but will
    occassionally and unexpectedly return an http 503 error when doing so.
    """
    from Bio import Entrez
    from Bio._py3k import HTTPError

    if isinstance(accessions, str):
        accessions = accessions.split(',')

    accessions = list(accessions)

    if email:
        Entrez.email = email

    head = 0
    tail = retmax = 10000

    try:
        while accessions[head:tail]:
            handle = Entrez.efetch(
                db='nucleotide',
                id=accessions[head:tail],
                retmode='xml',
                retmax=retmax,
                rettype='fasta')

            head = tail
            tail += retmax

            for record in Entrez.parse(handle):
                # split to remove accession version number
                yield (record['TSeq_accver'].split('.', 1)[0],
                       record['TSeq_taxid'])

    except HTTPError as err:
        # no need to crash everything if no results
        log.error(err)


def build_parser(parser):

    parser.add_argument(
        '--email',
        help='email address for use with Bio.Entrez.efetch')
    parser.add_argument(
        '--accession', metavar='ACCESSIONS',
        dest='accessions',
        help=('list of accession numbers provided as '
              'a comma-delimited list on the command line')
    )
    parser.add_argument(
        '--accession-file', metavar='FILE',
        dest='accessions_file', type=argparse.FileType('rU'),
        help=('file containing a list of accession numbers, one per line')
    )

    parser.add_argument(
        '-o', '--out-file', metavar='CSV', type=argparse.FileType('w'),
        dest='outfile',
        help='output file', default=sys.stdout
    )


def action(args):
    csv_out = csv.writer(args.outfile)

    accessions = set()
    if args.accessions:
        accessions |= set(args.accessions.split(','))

    if args.accessions_file:
        accessions |= set(a.rstrip() for a in args.accessions_file)

    if accessions:
        csv_out.writerows(id_from_accessions(accessions, email=args.email))
