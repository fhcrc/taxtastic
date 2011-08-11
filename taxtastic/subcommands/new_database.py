"""Creates a CSV file describing lineages for a set of taxa"""

from taxtastic import ncbi
import os
from os import path
import logging

log = logging.getLogger(__name__)

def build_parser(parser):

    parser.add_argument(
        '-d', '--database-file',
        dest = 'database_file',
        default = 'ncbi_taxonomy.db',
        metavar = 'FILE',
        help = """Name of the sqlite database file [%(default)s].""")

    parser.add_argument(
        '-p', '--download-dir',
        dest = 'download_dir',
        default = None,
        metavar = 'PATH',
        help = """Name of the directory into which to download the zip
        archive. [default is the same directory as the database file]""")

    parser.add_argument(
        '-x', '--clobber', action = 'store_true',
        dest = 'clobber', default = False,
        help = """Download a new zip archive containing NCBI taxonomy
        and/or re-create the database even if one or both already
        exists. [%(default)s]""")

def action(args):

    dbname = args.database_file
    pth, fname = path.split(dbname)
    zip_dest = args.download_dir or pth or '.'

    zfile, downloaded = ncbi.fetch_data(
        dest_dir = zip_dest,
        clobber = args.clobber)

    if not os.access(dbname, os.F_OK) or args.clobber:
        log.warning('creating new database in %s using data in %s' % \
                        (dbname, zfile))
        con = ncbi.db_connect(dbname, clobber=True)
        ncbi.db_load(con, zfile)
        con.close()
    else:
        log.warning('taxonomy database already exists in %s' % dbname)

