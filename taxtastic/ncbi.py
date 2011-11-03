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
"""
Methods and variables specific to the NCBI taxonomy.
"""

import sqlite3
import itertools
import logging
import zipfile

from . import taxdb
from .errors import IntegrityError
from .taxdb import do_insert, db_schema


log = logging

ncbi_data_url = 'ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdmp.zip'

# Schema specific to NCBI taxonomy

ranks = taxdb.ranks[:]

def db_connect(dbname='ncbi_taxonomy.db', schema=db_schema, clobber=False):
    """
    Create a connection object to a database. Attempt to establish a
    schema. If there are existing tables, delete them if clobber is
    True and return otherwise. Returns a connection object.
    """
    return taxdb.db_connect(dbname=dbname, schema=schema, clobber=clobber)

def db_load(con, archive, root_name='root', maxrows=None):
    """
    Load data from zip archive into database identified by con. Data
    is not loaded if target tables already contain data.
    """

    try:
        # nodes
        rows = read_nodes(
            rows=read_archive(archive, 'nodes.dmp'),
            root_name=root_name,
            ncbi_source_id=1)
        do_insert(con, 'nodes', rows, maxrows, add=False)

        # names
        rows = read_names(
            rows=read_archive(archive, 'names.dmp')
            )
        do_insert(con, 'names', rows, maxrows, add=False)

        # merged
        rows = read_archive(archive, 'merged.dmp')
        do_insert(con, 'merged', rows, maxrows, add=False)

    except sqlite3.IntegrityError, err:
        raise IntegrityError(err)

def fetch_data(dest_dir='.', clobber=False, url=ncbi_data_url):

    """
    Download data from NCBI required to generate local taxonomy
    database. Default url is ncbi.ncbi_data_url

    * dest_dir - directory in which to save output files (created if necessary).
    * clobber - don't download if False and target of url exists in dest_dir
    * url - url to archive; default is ncbi.ncbi_data_url

    Returns (fname, downloaded), where fname is the name of the
    downloaded zip archive, and downloaded is True if a new files was
    downloaded, false otherwise.

    see ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump_readme.txt
    """
    return taxdb.fetch_url(url, dest_dir, clobber)

def read_archive(archive, fname):
    """
    Return an iterator of rows from a zip archive.

    * archive - path to the zip archive.
    * fname - name of the compressed file within the archive.
    """

    zfile = zipfile.ZipFile(archive, 'r')
    for line in zfile.read(fname).splitlines():
        yield line.rstrip('\t|\n').split('\t|\t')

def read_dmp(fname):
    for line in open(fname,'rU'):
        yield line.rstrip('\t|\n').split('\t|\t')

def read_nodes(rows, root_name, ncbi_source_id):
    """
    Return an iterator of rows ready to insert into table "nodes".

    * rows - iterator of lists (eg, output from read_archive or read_dmp)
    * root_name - string identifying the root node (replaces NCBI's default).
    """

    keys = 'tax_id parent_id rank embl_code division_id'.split()
    idx = dict((k,i) for i,k in enumerate(keys))
    tax_id, parent_id, rank = [idx[k] for k in ['tax_id','parent_id','rank']]

    # assume the first row is the root
    row = rows.next()
    row[rank] = root_name
    rows = itertools.chain([row], rows)

    ncol = len(keys)
    # replace whitespace in "rank" with underscore
    for row in rows:
        row[rank] = '_'.join(row[rank].split())
        yield row[:ncol] + [ncbi_source_id]

def read_names(rows):
    """
    Return an iterator of rows ready to insert into table
    "names". Adds column "is_primary".

    * rows - iterator of lists (eg, output from read_archive or read_dmp)
    """

    keys = 'tax_id tax_name unique_name name_class'.split()
    idx = dict((k,i) for i,k in enumerate(keys))
    tax_name, unique_name, name_class = \
        [idx[k] for k in ['tax_name', 'unique_name', 'name_class']]

    def _is_primary(row):
        """
        Defines a name as "primary," meaning that other names associated
        with the this tax_id should be considered synonyms.
        """

        if row[name_class] != 'scientific name':
            result = 0
        elif not row[unique_name]:
            result = 1
        elif row[tax_name] == row[unique_name].split('<')[0].strip():
            result = 1
        else:
            result = 0

        return result

    # appends additional field is_primary
    for row in rows:
        yield row + [_is_primary(row)]


