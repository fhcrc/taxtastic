"""Find the intersection of a taxtable and a refpkg's taxonomy."""

import sqlite3
import logging
import sys
import csv
import argparse
import tempfile

from taxtastic.refpkg import Refpkg
from taxtastic.taxdb import Taxdb

log = logging.getLogger(__name__)

def build_parser(parser):
    parser.add_argument('infile', type=argparse.FileType('r'),
                        help='taxtable to compare against')
    parser.add_argument('-c', '--refpkg',
                        help='refpkg to insert into')
    parser.add_argument('-o','--outfile', type=argparse.FileType('w'), default=sys.stdout,
                        help='output file in csv format (default is stdout)')

def action(args):
    rp = Refpkg(args.refpkg)
    rp.load_db()
    cursor = rp.db.cursor()

    with tempfile.NamedTemporaryFile() as tmp_db:
        taxtable_db = Taxdb(sqlite3.connect(tmp_db.name))
        taxtable_db.create_tables()
        reader = csv.DictReader(args.infile)
        taxtable_db.insert_from_taxtable(lambda: reader._fieldnames, reader)
        cursor.execute('ATTACH DATABASE ? AS tt', (tmp_db.name,))

        writer = csv.writer(args.outfile)
        writer.writerow(('tax_id', 'intersection_rank'))
        cursor.execute("""
            SELECT child,
                   rank
            FROM   (SELECT child,
                           rank
                    FROM   parents
                           JOIN tt.taxa
                             ON tax_id = parent
                           JOIN tt.ranks USING (rank)
                    ORDER  BY child,
                              rank_order ASC)
            GROUP  BY child
        """)
        writer.writerows(cursor)
