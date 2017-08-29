"""Find the intersection of a taxtable and a refpkg's taxonomy."""

import logging
import sys
import csv
import argparse
import tempfile
import itertools
import operator

try:
    from pysqlite2 import dbapi2 as sqlite3
except ImportError:
    import sqlite3

from taxtastic.refpkg import Refpkg
from taxtastic.taxdb import Taxdb

log = logging.getLogger(__name__)


def build_parser(parser):
    parser.add_argument('infile', type=argparse.FileType('r'),
                        help='taxtable to compare against')
    parser.add_argument('-c', '--refpkg', required=True,
                        help='refpkg to insert into')
    parser.add_argument('-r', '--ranks', required=True,
                        help='ranks to list in the output')
    parser.add_argument('--all-ranks', default=False, action='store_true',
                        help="don't filter by the lowest rank; list all intersections")
    parser.add_argument('-o', '--out', type=argparse.FileType('w'), default=sys.stdout,
                        help='output file in csv format (default is stdout)')


def test_output(infile, outfile, ranks):
    """
    Ensure that all input taxids (provided in a taxtable) are represented in the
    output.
    """

    with open(infile, 'rU') as i, open(outfile, 'rU') as o:
        taxids_in = set(d['tax_id']
                        for d in csv.DictReader(i) if d['rank'] in ranks)
        taxids_out = set(d['tax_id'] for d in csv.DictReader(o))

        assert len(taxids_in - taxids_out) == 0, taxids_in - taxids_out


def filter_ranks(results):
    """
    Find just the first rank for all the results for a given tax_id.
    """
    for _, group in itertools.groupby(results, operator.itemgetter(0)):
        yield next(group)


def action(args):
    rp = Refpkg(args.refpkg, create=False)
    rp.load_db()
    cursor = rp.db.cursor()
    ranks = args.ranks.split(',')

    with tempfile.NamedTemporaryFile() as tmp_db:
        taxtable_db = Taxdb(sqlite3.connect(tmp_db.name))
        taxtable_db.create_tables()
        reader = csv.DictReader(args.infile)
        taxtable_db.insert_from_taxtable(lambda: reader._fieldnames, reader)
        cursor.execute('ATTACH DATABASE ? AS tt', (tmp_db.name,))

        writer = csv.writer(args.out)
        writer.writerow(('tax_id', 'intersection_rank'))
        cursor.execute("""
            SELECT tax_id,
                   COALESCE(itaxa.rank, "")
              FROM tt.taxa
                   LEFT JOIN (SELECT child AS tax_id,
                                     rank_order,
                                     rank
                                FROM tt.parents
                                     JOIN taxa
                                       ON tax_id = parent
                                     JOIN ranks USING (rank)
                               WHERE rank IN (%s)) itaxa USING (tax_id)
             ORDER BY tax_id,
                      rank_order DESC
        """ % ', '.join('?' * len(ranks)), ranks)
        if not args.all_ranks:
            cursor = filter_ranks(cursor)
        writer.writerows(cursor)

    args.out.flush()
    test_output(args.infile.name, args.out.name, ranks)
