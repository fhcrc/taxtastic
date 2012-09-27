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
import os
import urllib
import zipfile
import re

from errors import IntegrityError

log = logging

ncbi_data_url = 'ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdmp.zip'

db_schema = """
-- nodes.dmp specifies additional columns but these are not implemented yet
CREATE TABLE nodes(
tax_id        TEXT UNIQUE PRIMARY KEY NOT NULL,
parent_id     TEXT,
rank          TEXT,
embl_code     TEXT,
division_id   INTEGER,
source_id     INTEGER DEFAULT 1 -- added to support multiple sources
);

CREATE TABLE names(
tax_id        TEXT REFERENCES nodes(tax_id),
tax_name      TEXT,
unique_name   TEXT,
name_class    TEXT,
-- not defined in names.dmp:
is_primary    INTEGER,
is_classified INTEGER -- tax_name does not match UNCLASSIFIED_REGEX
);

CREATE TABLE merged(
old_tax_id    TEXT,
new_tax_id    TEXT REFERENCES nodes(tax_id)
);

-- table "source" supports addition of custom taxa (not provided by NCBI)
CREATE TABLE source(
id            INTEGER PRIMARY KEY AUTOINCREMENT,
name          TEXT UNIQUE,
description   TEXT
);

INSERT INTO "source"
  (id, name, description)
VALUES
  (1, "NCBI", "NCBI taxonomy");

-- indices on nodes
CREATE INDEX nodes_tax_id ON nodes(tax_id);
CREATE INDEX nodes_parent_id ON nodes(parent_id);
CREATE INDEX nodes_rank ON nodes(rank);

-- indices on names
CREATE INDEX names_tax_id ON names(tax_id);
CREATE INDEX names_tax_name ON names(tax_name);
CREATE INDEX names_is_primary ON names(is_primary);
CREATE INDEX names_is_classified ON names(is_classified);
CREATE INDEX names_taxid_is_primary ON names(tax_id, is_primary);
CREATE INDEX names_name_is_primary ON names(tax_name, is_primary);
-- CREATE UNIQUE INDEX names_id_name ON names(tax_id, tax_name, is_primary);

"""

# define headers in names.dmp, etc (may not correspond to table columns above)
merged_keys = 'old_tax_id new_tax_id'.split()

undefined_rank = 'no_rank'
root_name = 'root'

_ranks = """
root
superkingdom
kingdom
subkingdom
superphylum
phylum
subphylum
superclass
class
subclass
infraclass
superorder
order
suborder
infraorder
parvorder
superfamily
family
subfamily
tribe
subtribe
genus
subgenus
species group
species subgroup
species
subspecies
varietas
forma
"""

# provides criteria for defining matching tax_ids as "unclassified"
UNCLASSIFIED_REGEX = re.compile(
    r'' + r'|'.join(frozenset(['-like',
                               'Taxon'
                               '\d\d',
                               'acidophile',
                               'actinobacterium',
                               'aerobic',
                               r'\b[Al]g(um|a)\b',
                               r'\b[Bb]acteri(um|a)',
                               'Barophile',
                               'cyanobacterium',
                               'Chloroplast',
                               'Cloning',
                               'cluster',
                               '-containing',
                               'epibiont',
                               # 'et al',
                               'eubacterium',
                               r'\b[Gg]roup\b',
                               'halophilic',
                               r'hydrothermal\b',
                               'isolate',
                               'marine',
                               'methanotroph',
                               'microorganism',
                               'mollicute',
                               'pathogen',
                               '[Pp]hytoplasma',
                               'proteobacterium',
                               'putative',
                               r'\bsp\.',
                               'species',
                               'spirochete',
                               r'str\.'
                               'strain',
                               'symbiont',
                               'taxon',
                               'unicellular',
                               'uncultured',
                               'unclassified',
                               'unidentified',
                               'unknown',
                               'vector\b',
                               r'vent\b',
                               ])))

ranks = [k.strip().replace(' ','_') for k in _ranks.splitlines() if k.strip()]

def db_connect(dbname='ncbi_taxonomy.db', schema=db_schema, clobber = False):
    """
    Create a connection object to a database. Attempt to establish a
    schema. If there are existing tables, delete them if clobber is
    True and return otherwise. Returns a connection object.
    """

    if clobber:
        log.info('Creating new database %s' % dbname)
        try:
            os.remove(dbname)
        except OSError:
            pass

    con = sqlite3.connect(dbname)
    cur = con.cursor()

    cmds = [cmd.strip() for cmd in schema.split(';') if cmd.strip()]
    try:
        for cmd in cmds:
            cur.execute(cmd)
            log.debug(cmd)
    except sqlite3.OperationalError as err:
        log.info(err)

    return con

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
            rows=read_archive(archive, 'names.dmp'),
            unclassified_regex = UNCLASSIFIED_REGEX
            )
        do_insert(con, 'names', rows, maxrows, add=False)

        # merged
        rows = read_archive(archive, 'merged.dmp')
        do_insert(con, 'merged', rows, maxrows, add=False)

        fix_missing_primary(con)

    except sqlite3.IntegrityError, err:
        raise IntegrityError(err)

def fix_missing_primary(con):
    missing_primary = """SELECT tax_id
        FROM names
        GROUP BY tax_id
        HAVING SUM(is_primary) = 0;"""
    rows_for_taxid = """SELECT tax_id, tax_name, unique_name, name_class
        FROM names
        WHERE tax_id = ?"""
    cursor = con.cursor()
    tax_ids = [i[0] for i in cursor.execute(missing_primary)]
    logging.warn("%d records lack primary names", len(tax_ids))

    for tax_id in tax_ids:
        records = list(cursor.execute(rows_for_taxid, [tax_id]))
        # Prefer scientific name

        if sum(i[-1] == 'scientific name' for i in records) == 1:
            tax_id, tax_name, unique_name, name_class = next(i for i in records
                    if i[-1] == 'scientific name')
            cursor.execute("""UPDATE names
                SET is_primary = 1
                WHERE tax_id = ? AND name_class = ?""",
                [tax_id, 'scientific name'])
        else:
            tax_id, tax_name, unique_name, name_class = records[0]
        logging.warn("No primary name for tax_id %s. Arbitrarily using %s[%s].",
                tax_id, tax_name, name_class)
        cursor.execute("""UPDATE names
            SET is_primary = 1
            WHERE tax_id = ? AND tax_name = ? AND
                unique_name = ? AND name_class = ?""",
            [tax_id, tax_name, unique_name, name_class])

def do_insert(con, tablename, rows, maxrows=None, add=True):
    """
    Insert rows into a table. Do not perform the insert if
    add is False and table already contains data.
    """

    cur = con.cursor()

    cur.execute('select count(*) from "%s" where rowid < 2' % tablename)
    has_data = cur.fetchone()[0]

    if not add and has_data:
        log.info('Table "%s" already contains data; load not performed.' % tablename)
        return False

    # pop first row to determine number of columns
    row = rows.next()
    cmd = 'INSERT INTO "%s" VALUES (%s)' % (tablename, ', '.join(['?']*len(row)))
    log.info(cmd)

    # put the first row back
    rows = itertools.chain([row], rows)
    if maxrows:
        rows = itertools.islice(rows, maxrows)

    cur.executemany(cmd, rows)
    con.commit()

    return True

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

    dest_dir = os.path.abspath(dest_dir)
    try:
        os.mkdir(dest_dir)
    except OSError:
        pass

    fout = os.path.join(dest_dir, os.path.split(url)[-1])

    if os.access(fout, os.F_OK) and not clobber:
        downloaded = False
        log.warning('%s exists; not downloading' % fout)
    else:
        downloaded = True
        log.warning('downloading %(url)s to %(fout)s' % locals())
        urllib.urlretrieve(url, fout)

    return (fout, downloaded)

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

def read_names(rows, unclassified_regex = None):
    """
    Return an iterator of rows ready to insert into table
    "names". Adds columns "is_primary" and "is_classified". If
    `unclassified_regex` is not None, defines 'is_classified' as 1 if
    the regex fails to match "tax_name" or 0 otherwise; if
    `unclassified_regex` is None, 'is_classified' is given a value of
    None.

    * rows - iterator of lists (eg, output from read_archive or read_dmp)
    * unclassified_regex - a compiled re matching "unclassified" names
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

    if unclassified_regex:
        def _is_classified(row):
            """
            Return 1 if tax_name element of `row` matches
            unclassified_regex, 0 otherwise. Search no more than the
            first two whitespace-delimited words.
            """
            tn = row[tax_name]
            term = ' '.join(tn.split()[:2]) if ' ' in tn else tn
            return 0 if unclassified_regex.search(term) else 1
    else:
        _is_classified = lambda row: None

    # appends additional field is_primary
    for row in rows:
        yield row + [_is_primary(row), _is_classified(row)]
