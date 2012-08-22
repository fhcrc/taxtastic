"""Create a database containing an entire taxonomy"""
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

    parser.add_argument(
        '--preserve-inconsistent-taxonomies',
        action='store_true', default=False,
        help="""If a node has the same rank as its parent, do *not* its rank
        set to no_rank.""")

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
        with con:
            ncbi.db_load(con, zfile)
            if not args.preserve_inconsistent_taxonomies:
                curs = con.cursor()
                curs.execute("""
                    UPDATE nodes
                       SET rank = 'no_rank'
                     WHERE tax_id IN (SELECT n1.tax_id
                                        FROM nodes n1
                                             JOIN nodes n2
                                               ON n1.parent_id = n2.tax_id
                                       WHERE n1.rank = n2.rank
                                         AND n1.rank NOT IN ('root', 'no_rank'))
                """)
        con.close()
    else:
        log.warning('taxonomy database already exists in %s' % dbname)

