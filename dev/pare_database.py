#!/usr/bin/env python
"""
Generate a subset of the NCBI taxonomy database to facilitate testing

Generate a new copy with:
    rm -f testfiles/small_taxonomy.db
    python devtools/pare_database.py -k testfiles/keep_taxids.txt | sqlite3 testfiles/small_taxonomy.db
"""

import argparse
import re
import sqlite3
import shutil
import sys
import tempfile

import sqlalchemy

from taxtastic import taxonomy, ncbi


def in_clause(args):
    return '({0})'.format(', '.join('?' * len(args)))


def lineage(taxonomy, tax_id):
    try:
        return taxonomy._get_lineage(tax_id=tax_id)
    except KeyError:
        m = taxonomy._get_merged(tax_id)
        if m != tax_id:
            return lineage(taxonomy, m)
        else:
            print >> sys.stderr, "unknown taxid:", tax_id
            return []


def main():
    parser = argparse.ArgumentParser(description="""Generate SQL statements to
            create a smaller version of a taxonomy database, keeping only
            lineages associated with a set of tax_ids Use as:

            python devtools/pare_database.py -k testfiles/keep_taxids.txt | sqlite3 test_output/ncbi_master.db""")

    parser.add_argument('-d', help='source database',
                        default='ncbi_taxonomy.db')
    parser.add_argument('-k', '--keep-taxids', help="""file containing
        whitespace-delimited list of taxids to keep""",
                        type=argparse.FileType('r'))
    parser.add_argument('-o', '--output-file',
                        default=sys.stdout, type=argparse.FileType('w'))

    a = parser.parse_args()
    with a.keep_taxids as fp:
        lines = (i for i in fp if not i.startswith('#'))
        keep_ids = frozenset(i for i in re.split(
            r'[\s\n]+', ''.join(lines), flags=re.MULTILINE) if i)

    e = sqlalchemy.create_engine('sqlite:///{0}'.format(a.d))
    t = taxonomy.Taxonomy(e, ncbi.ranks)

    # Get lineages
    lineages = (lineage(t, i) for i in keep_ids)
    keep_taxa = frozenset(i[1] for l in lineages for i in l)

    with tempfile.NamedTemporaryFile() as tf:
        with open(a.d) as fp:
            print >> sys.stderr, "copying db"
            shutil.copyfileobj(fp, tf, 20 << 10)
        tf.flush()

        con = sqlite3.connect(tf.name)
        cur = con.cursor()
        ic = in_clause(keep_taxa)
        print >> sys.stderr, "pruning nodes"
        cur.execute("DELETE FROM nodes WHERE tax_id NOT IN {0}".format(
            ic), list(keep_taxa))
        print >> sys.stderr, "pruning names"
        cur.execute("DELETE FROM names WHERE tax_id NOT IN {0}".format(
            ic), list(keep_taxa))
        print >> sys.stderr, "pruning merged"
        cur.execute("DELETE FROM merged WHERE old_tax_id NOT IN {0} AND new_tax_id NOT IN {0}".format(ic),
                    list(keep_taxa) + list(keep_taxa))

        with a.output_file as fp:
            for stmt in con.iterdump():
                print >> fp, stmt


if __name__ == '__main__':
    main()
