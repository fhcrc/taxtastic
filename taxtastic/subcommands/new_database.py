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
already exists it be will overwritten unless you specify ``--no-clobber``.
The NCBI taxonomy will be downloaded into
the same directory as ``database_file`` will be created in unless you
specify ``-p`` or ``--download-dir``.
"""
import argparse
import logging
import sqlalchemy
import sys

import taxtastic

log = logging.getLogger(__name__)


def build_parser(parser):
    parser = taxtastic.utils.add_database_args(parser)

    parser.add_argument(
        '--no-clobber',
        action='store_false', default=True,
        dest='clobber',
        help=('If database exists keep current data '
              'and append new data. [False]'))

    download_parser = parser.add_argument_group(title='download options')
    download_parser.add_argument(
        '-z', '--taxdump-file',
        metavar='FILE.zip',
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

    parser.add_argument(
        '--out',
        type=argparse.FileType('w'),
        default=sys.stdout,
        help='table sql')


def action(args):
    if args.taxdump_file:
        zfile = args.taxdump_file
    else:
        zip_dest = args.download_dir or '.'
        zfile, _ = taxtastic.ncbi.fetch_data(
            dest_dir=zip_dest,
            clobber=args.clobber,
            url=args.taxdump_url)

    engine = sqlalchemy.create_engine(args.url, echo=args.verbosity > 2)

    # creates database schema
    taxtastic.ncbi.db_connect(engine, schema=args.schema, clobber=args.clobber)

    ncbi_loader = taxtastic.ncbi.NCBILoader(engine, args.schema)
    ncbi_loader.load_archive(zfile)
    ncbi_loader.set_names_is_classified()
    ncbi_loader.set_nodes_is_valid()

    # print_sql(args.out, engine.name, base.metadata)


def print_sql(out, engine_name, metadata):
    def dump(sql, *multiparams, **params):
        out.write(str(sql.compile(dialect=dump.dialect)).strip() + ';\n')
    engine = sqlalchemy.create_engine(
        engine_name + '://', strategy='mock', executor=dump)

    dump.dialect = engine.dialect
    metadata.create_all(bind=engine)
