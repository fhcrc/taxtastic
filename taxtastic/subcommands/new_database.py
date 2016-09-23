"""Download NCBI taxonomy and create a database

Download the current version of the NCBI taxonomy and load it into
``database_file`` as an SQLite3 database.  If ``database_file``
already exists, it will fail and leave it untouched unless you specify
``-x`` or ``--clobber``.  The NCBI taxonomy will be downloaded into
the same directory as ``database_file`` will be created in unless you
specify ``-p`` or ``--download-dir``.

"""
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

import os
import logging
import taxtastic

log = logging.getLogger(__name__)


def build_parser(parser):

    parser.add_argument(
        '-d', '--database-file',
        dest='database_file',
        default='ncbi_taxonomy.db',
        metavar='FILE',
        help="""Name of the sqlite database file [%(default)s].""")

    parser.add_argument(
        '-z', '--taxdump-file',
        metavar='ZIP',
        help='Location of zipped taxdump file [%(default)s]')

    parser.add_argument(
        '-u', '--taxdump-url',
        default=taxtastic.ncbi.DATA_URL,
        metavar='URL',
        help='Url to taxdump file [%(default)s]')

    parser.add_argument(
        '-p', '--download-dir',
        dest='download_dir',
        metavar='PATH',
        help="""Name of the directory into which to download the zip
             archive. [default is the same directory as the database file]""")

    parser.add_argument(
        '-x', '--clobber', action='store_true',
        dest='clobber', default=False,
        help="""Download a new zip archive containing NCBI taxonomy
        and/or re-create the database even if one or both already
        exists. [%(default)s]""")

    parser.add_argument(
        '--preserve-inconsistent-taxonomies',
        action='store_true', default=False,
        help="""If a node has the same rank as its parent, do *not* its rank
        set to no_rank. [%(default)s]""")


def action(args):
    dbname = args.database_file

    if not os.access(args.database_file, os.F_OK) or args.clobber:
        if args.taxdump_file:
            zfile = args.taxdump_file
        else:
            pth = os.path.split(dbname)[0]
            zip_dest = args.download_dir or pth or '.'
            zfile, _ = taxtastic.ncbi.fetch_data(
                dest_dir=zip_dest,
                clobber=args.clobber,
                url=args.taxdump_url)
        msg = 'creating new database in {} using data in {}'
        log.warning(msg.format(dbname, zfile))
        engine = taxtastic.ncbi.db_connect(dbname, clobber=args.clobber)
        taxtastic.ncbi.db_load(engine, zfile)
    else:
        log.warning('taxonomy database already exists in %s' % dbname)
