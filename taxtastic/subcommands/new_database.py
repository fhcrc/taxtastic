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
"""Download NCBI taxonomy and create a database

Download the current version of the NCBI taxonomy and load it into
``database_file`` as an SQLite3 database.  If ``database_file``
already exists, it will fail and leave it untouched unless you specify
``-x`` or ``--clobber``.  The NCBI taxonomy will be downloaded into
the same directory as ``database_file`` will be created in unless you
specify ``-p`` or ``--download-dir``.
"""
import logging
import taxtastic

log = logging.getLogger(__name__)


def build_parser(parser):

    parser.add_argument(
        'url',
        default='sqlite:///ncbi_taxonomy.db',
        help='url to database [%(default)s]')
    parser.add_argument(
        '--schema',
        help='database schema to use if applicable')
    parser.add_argument(
        '--append',
        action='store_false',
        dest='clobber',
        help=('If database exists keep current data '
              'and append new data. [False]'))

    download_parser = parser.add_argument_group(title='download options')
    download_parser.add_argument(
        '-z', '--taxdump-file',
        metavar='ZIP',
        help='Location of zipped taxdump file [taxdmp.zip]')

    download_parser.add_argument(
        '-u', '--taxdump-url',
        default=taxtastic.ncbi.DATA_URL,
        metavar='URL',
        help='Url to taxdump file [%(default)s]')

    download_parser.add_argument(
        '-p', '--download-dir',
        dest='download_dir',
        metavar='PATH',
        help="""Name of the directory into which to download the zip
             archive. [default is the same directory as the database file]""")


def action(args):
    if args.taxdump_file:
        zfile = args.taxdump_file
    else:
        zip_dest = args.download_dir or '.'
        zfile, _ = taxtastic.ncbi.fetch_data(
            dest_dir=zip_dest,
            clobber=args.clobber,
            url=args.taxdump_url)
    engine = taxtastic.ncbi.db_connect(
        args.url,
        schema=args.schema,
        clobber=args.clobber,
        verbosity=args.verbosity)
    taxtastic.ncbi.db_load(engine, zfile, schema=args.schema)
