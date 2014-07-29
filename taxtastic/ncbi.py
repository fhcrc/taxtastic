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

import itertools
import logging
import operator
import os
import re
import sqlite3
import urllib
import zipfile

import sqlalchemy
from sqlalchemy import Column, Integer, String, Boolean, ForeignKey, Index
from sqlalchemy.orm import relationship
from sqlalchemy.ext.declarative import declarative_base

from errors import IntegrityError

log = logging

Base = declarative_base()

ncbi_data_url = 'ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdmp.zip'

class Node(Base):
    __tablename__ = 'nodes'
    tax_id = Column(String, primary_key=True, nullable=False)
    parent_id = Column(String, index=True)
    rank = Column(String, index=True)
    embl_code = Column(String)
    division_id = Column(String)
    source_id = Column(Integer, server_default='1')
    is_valid = Column(Boolean, server_default='1', index=True)

class Name(Base):
    __tablename__ = 'names'
    id = Column(Integer, primary_key=True)

    tax_id = Column(String, ForeignKey('nodes.tax_id', ondelete='CASCADE'), index=True)
    node = relationship('Node', backref='names')
    tax_name = Column(String, index=True)
    unique_name = Column(String)
    name_class = Column(String)
    is_primary = Column(Boolean)
    is_classified = Column(Boolean)
Index('ix_names_tax_id_is_primary', Name.tax_id, Name.is_primary)

class Merge(Base):
    __tablename__ = 'merged'

    old_tax_id = Column(String, primary_key=True, index=True)
    new_tax_id = Column(String, ForeignKey('nodes.tax_id', ondelete='CASCADE'), index=True)
    merged_node = relationship('Node', backref='merged_ids')

class Source(Base):
    __tablename__ = 'source'

    id = Column(Integer, primary_key=True)
    name = Column(String, unique=True)
    description = Column(String)


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

# Components of a regex to apply to all names. Names matching this regex are
# marked as invalid.
UNCLASSIFIED_REGEX_COMPONENTS = [r'-like\b',
                                r'\bactinomycete\b',
                                r'\bcrenarchaeote\b',
                                r'\bculture\b',
                                r'\bchimeric\b',
                                r'\bcyanobiont\b',
                                r'degrading',
                                r'\beuryarchaeote\b',
                                r'disease',
                                r'\b[cC]lone',
                                r'\bmethanogen(ic)?\b',
                                r'\bplanktonic\b',
                                r'\bplanctomycete\b',
                                r'\bsymbiote\b',
                                r'\btransconjugant\b',
                                r'^[a-z]', # starts with lower-case character
                                r'^\W+\s+[a-zA-Z]*\d', # Digit in second word
                                r'\d\d',
                                r'atypical',
                                r'^cf\.',
                                r'acidophile',
                                r'\bactinobacterium\b',
                                r'aerobic',
                                r'.+\b[Al]g(um|a)\b',
                                r'\b[Bb]acteri(um|al)\b',
                                r'.+\b[Bb]acteria\b',
                                r'Barophile',
                                r'cyanobacterium',
                                r'Chloroplast',
                                r'Cloning',
                                r'\bclone\b',
                                r'cluster',
                                r'^diazotroph',
                                r'\bcoccus\b',
                                r'archaeon',
                                r'-containing',
                                r'epibiont',
                                # 'et al',
                                r'environmental samples',
                                r'eubacterium',
                                r'\b[Gg]roup\b',
                                r'halophilic',
                                r'hydrothermal\b',
                                r'isolate',
                                r'\bmarine\b',
                                r'methanotroph',
                                r'microorganism',
                                r'mollicute',
                                r'pathogen',
                                r'[Pp]hytoplasma',
                                r'proteobacterium',
                                r'putative',
                                r'\bsp\.',
                                r'species',
                                r'spirochete',
                                r'str\.',
                                r'strain',
                                r'symbiont',
                                r'\b[Tt]axon\b',
                                r'unicellular',
                                r'uncultured',
                                r'unclassified',
                                r'unidentified',
                                r'unknown',
                                r'vector\b',
                                r'vent\b',
                               ]

# provides criteria for defining matching tax_ids as "unclassified"
UNCLASSIFIED_REGEX = re.compile('|'.join(UNCLASSIFIED_REGEX_COMPONENTS))

ranks = [k.strip().replace(' ','_') for k in _ranks.splitlines() if k.strip()]

def db_connect(dbname='ncbi_taxonomy.db', clobber=False):
    """
    Create a connection object to a database. Attempt to establish a
    schema. If there are existing tables, delete them if clobber is
    True and return otherwise. Returns a sqlalchemy engine object.
    """

    if clobber:
        log.info('Creating new database %s' % dbname)
        try:
            os.remove(dbname)
        except OSError:
            pass

    engine = sqlalchemy.create_engine('sqlite:///{0}'.format(dbname))
    Base.metadata.create_all(bind=engine)
    return engine

def db_load(engine, archive, root_name='root', maxrows=None):
    """
    Load data from zip archive into database identified by con. Data
    is not loaded if target tables already contain data.
    """

    try:
        # nodes
        logging.info("Inserting nodes")
        rows = read_nodes(
            rows=read_archive(archive, 'nodes.dmp'),
            root_name=root_name,
            ncbi_source_id=1)
        # Add is_valid
        do_insert(engine, 'nodes', rows, maxrows, add=False)

        # names
        logging.info("Inserting names")
        rows = read_names(
            rows=read_archive(archive, 'names.dmp'),
            unclassified_regex=UNCLASSIFIED_REGEX)
        do_insert(engine, 'names', rows, maxrows, add=False)

        # merged
        logging.info("Inserting merged")
        rows = read_archive(archive, 'merged.dmp')
        rows = (dict(zip(['old_tax_id', 'new_tax_id'], row)) for row in rows)
        do_insert(engine, 'merged', rows, maxrows, add=False)

        fix_missing_primary(engine)

        # Mark names as valid/invalid
        mark_is_valid(engine)
        update_subtree_validity(engine)

    except sqlite3.IntegrityError, err:
        raise IntegrityError(err)

def fix_missing_primary(engine):
    with engine.begin() as cursor:
        missing_primary = """SELECT tax_id
            FROM names
            GROUP BY tax_id
            HAVING SUM(is_primary) = 0;"""
        rows_for_taxid = """SELECT tax_id, tax_name, unique_name, name_class
            FROM names
            WHERE tax_id = ?"""
        tax_ids = [i[0] for i in cursor.execute(missing_primary)]
        logging.warn("%d records lack primary names", len(tax_ids))

        for tax_id in tax_ids:
            records = list(cursor.execute(rows_for_taxid, [tax_id]))
            # Prefer scientific name
            if sum(list(i)[-1] == 'scientific name' for i in records) == 1:
                tax_id, tax_name, unique_name, name_class = next(i for i in records
                        if list(i)[-1] == 'scientific name')
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

def mark_is_valid(engine, regex=UNCLASSIFIED_REGEX):
    """
    Apply ``regex`` to primary names associated with tax_ids, marking those
    that match as invalid.
    """
    logging.info("Marking nodes validity based on primary name")
    with engine.begin() as connection:
        sql = """UPDATE nodes SET is_valid = (SELECT is_classified FROM names WHERE names.tax_id = nodes.tax_id AND names.is_primary = 1)"""
        connection.execute(sql)

def partition(iterable, size):
    iterable = iter(iterable)
    while True:
        chunk = list(itertools.islice(iterable, 0, size))
        if chunk:
            yield chunk
        else:
            break

def update_subtree_validity(engine, mark_below_rank='species'):
    """
    Update subtrees below rank "species" to match ``is_valid`` status at
    rank "species"

    Also covers the special case of marking the "unclassified Bacteria" subtree
    invalid.
    """
    def generate_in_param(count):
        return '(' + ', '.join('?' * count) + ')'

    def mark_subtrees(conn, tax_ids, is_valid):
        to_mark = list(tax_ids)
        logging.info("Marking %d subtrees as is_valid=%s", len(to_mark), is_valid)
        while to_mark:
            # First, mark nodes
            conn.execute("""UPDATE nodes SET is_valid = ?
                WHERE tax_id = ?""",
                [[is_valid, tax_id] for tax_id in to_mark])

            # Find children - can exceed the maximum number of parameters in a sqlite query,
            # so chunk:
            chunked = partition(to_mark, 250)
            child_sql = """SELECT tax_id
                           FROM nodes
                           WHERE parent_id IN {0}"""
            to_mark = [i[0] for i in itertools.chain.from_iterable(
                           conn.execute(child_sql.format(generate_in_param(len(chunk))),
                                chunk) for chunk in chunked)]

    below_rank_query = """
    SELECT nodes.tax_id, pnodes.is_valid
    FROM nodes
        JOIN nodes pnodes ON pnodes.tax_id = nodes.parent_id
    WHERE pnodes.rank = ?
    ORDER BY pnodes.is_valid"""

    with engine.begin() as conn:
        result = list(conn.execute(below_rank_query, [mark_below_rank]))

        # Group by validity
        grouped = itertools.groupby(result, operator.itemgetter(1))
        for is_valid, records in grouped:
            tax_ids = [i[0] for i in records]
            mark_subtrees(conn, tax_ids, is_valid)

        # Special case: unclassified bacteria
        result = list(conn.execute("""SELECT tax_id FROM names WHERE tax_name = ? and is_primary = ?""",
            ['unclassified Bacteria', 1]))
        assert len(result) < 2
        logging.info("marking subtrees for unclassified Bacteria invalid")
        for i, in result:
            mark_subtrees(conn, [i], 0)


def do_insert(engine, tablename, rows, maxrows=None, add=True, chunk_size=5000):
    """
    Insert rows into a table. Do not perform the insert if
    add is False and table already contains data.
    """
    meta = Base.metadata
    table = meta.tables[tablename]
    has_data = table.select(bind=engine).limit(1).count().execute().first()[0] > 0

    if not add and has_data:
        log.info('Table "%s" already contains data; load not performed.' % tablename)
        return False
    if maxrows:
        rows = itertools.islice(rows, maxrows)

    insert = table.insert()
    with engine.begin() as conn:
        count = 0
        for chunk in partition(rows, chunk_size):
            result = conn.execute(insert, chunk)
            count += result.rowcount
        logging.info("Inserted %d rows into %s", count, tablename)

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
    for line in open(fname, 'rU'):
        yield line.rstrip('\t|\n').split('\t|\t')

def read_nodes(rows, root_name, ncbi_source_id):
    """
    Return an iterator of rows ready to insert into table "nodes".

    * rows - iterator of lists (eg, output from read_archive or read_dmp)
    * root_name - string identifying the root node (replaces NCBI's default).
    """

    keys = 'tax_id parent_id rank embl_code division_id'.split()
    idx = dict((k,i) for i,k in enumerate(keys))
    rank = idx['rank']

    # assume the first row is the root
    row = rows.next()
    row[rank] = root_name
    rows = itertools.chain([row], rows)

    colnames = keys + ['source_id']
    for r in rows:
        row = dict(zip(colnames, r))
        assert len(row) == len(colnames)

        # replace whitespace in "rank" with underscore
        row['rank'] = '_'.join(row['rank'].split())
        yield row

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
            return 0 if unclassified_regex.search(tn) else 1
    else:
        _is_classified = lambda row: None

    # appends additional field is_primary
    colnames = keys + ['is_primary', 'is_classified']
    for r in rows:
        row = dict(zip(colnames, r + [_is_primary(r), _is_classified(r)]))
        yield row
