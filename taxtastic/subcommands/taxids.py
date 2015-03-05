"""Look up a set of tax_ids from taxonomic names"""
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

from sqlalchemy import create_engine

from taxtastic.taxonomy import Taxonomy
from taxtastic.ncbi import ranks as ncbi_ranks

log = logging.getLogger(__name__)


def get_children(engine, parent_ids, rank='species'):
    """
    Recursively fetch children of tax_ids in `parent_ids` until the
    rank of `rank`
    """

    if not parent_ids:
        return []

    cmd = """
    select tax_id, tax_name, rank from nodes join names using (tax_id)
    where parent_id = "%s"
    and is_primary = 1
    """

    species = []
    for parent_id in parent_ids:
        result = engine.execute(cmd % parent_id)
        keys = result.keys()
        rows = [dict(zip(keys, row)) for row in result.fetchall()]
        species.extend(
            [r for r in rows if r['rank'] == rank
                and 'sp.' not in r['tax_name']])

        others = [r for r in rows if r['rank'] not in (rank, 'no_rank')]
        if others:
            _, s = get_children(engine, [r['tax_id'] for r in others])
            species.extend(s)

    return keys, species


def id_from_accessions(accessions, email=None):
    """
    Retrieve tax ids from a list of accessions.

    ncbi.efetch can work unexpectedly when querying a large number of tax ids.
    The best way to mitigate is to send it only its documented max of 10k ids
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
        '--database-file', dest='dbfile',
        action='store',
        help='Filename of sqlite database [%(default)s].',
        metavar='FILE')
    parser.add_argument(
        '-f', '--name-file', metavar='FILE',
        type=argparse.FileType('rU'),
        dest='taxnames_file',
        help=('file containing a list of taxonomic names, one per line')
    )
    parser.add_argument(
        '-n', '--name', metavar='NAMES',
        dest='taxnames',
        help=('list of taxonomic names provided as '
              'a comma-delimited list on the command line')
    )
    parser.add_argument(
        '--no-children', action='store_true',
        help='return only immediate tax_id ignoring species level children')

    output_group = parser.add_argument_group(
        "Output options").add_mutually_exclusive_group()
    output_group.add_argument(
        '-o', '--out-file', metavar='CSV', type=argparse.FileType('w'),
        dest='outfile',
        help='output file', default=sys.stdout
    )


def action(args):
    csv_out = csv.writer(args.outfile)

    if args.dbfile:
        dbfile = args.dbfile
        taxnames_file = args.taxnames_file
        taxnames = args.taxnames

        engine = create_engine('sqlite:///{}'.format(dbfile), echo=False)
        tax = Taxonomy(engine, ncbi_ranks)

        names = []
        if taxnames_file:
            names += [line.split('#', 1)[0].strip()
                      for line in taxnames_file if line.strip()
                      and not line.startswith('#')]

        if taxnames:
            names += [x.strip() for x in taxnames.split(',')]

        for name in set(names):
            tax_id, tax_name, is_primary, rank, note = '', '', '', '', ''

            try:
                tax_id, tax_name, is_primary = tax.primary_from_name(name)
            except KeyError:
                note = 'not found'
            else:
                parent, rank = tax._node(tax_id)
                note = '' if is_primary else 'not primary'

            if note:
                msg = '{name:>20} | {tax_id:>7} {tax_name:>20} {note}'
                msg = msg.format(**locals())
                log.warn(msg)

            if args.no_children or rank == 'species':
                csv_out.writerow((name, tax_name, tax_id))
            else:
                _, rows = get_children(engine, [tax_id], rank='species')
                for r in rows:
                    csv_out.writerow((name, r['tax_name'], r['tax_id']))
    elif args.taxnames or args.taxnames_file:
        log.error('no taxonomy database file specified')

    accessions = set()
    if args.accessions:
        accessions |= set(args.accessions.split(','))

    if args.accessions_file:
        accessions |= set(a.rstrip() for a in args.accessions_file)

    if accessions:
        csv_out.writerows(id_from_accessions(accessions, email=args.email))
