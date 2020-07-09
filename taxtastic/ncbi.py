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
import os
import re
from six.moves.urllib import request
# import urllib.parse
# import urllib.error
import zipfile
import io
from operator import itemgetter

import sqlalchemy
from sqlalchemy import (Column, Integer, String, Boolean,
                        ForeignKey, Index, MetaData, PrimaryKeyConstraint)
from sqlalchemy.orm import relationship
from sqlalchemy.ext.declarative import declarative_base

from taxtastic.utils import random_name


log = logging.getLogger(__name__)


DATA_URL = 'ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdmp.zip'

# For rank order: https://en.wikipedia.org/wiki/Taxonomic_rank
RANK_ORDER = [
    'forma',
    'subvariety',
    'varietas',
    'serogroup',
    'pathogroup',
    'morph',
    'genotype',
    'biotype',
    'subspecies',
    'species',
    'species_subgroup',
    'species_group',
    'subseries',
    'series',
    'subsection',
    'section',
    'subgenus',
    'genus',
    'subtribe',
    'tribe',
    'subfamily',
    'family',
    'superfamily',
    'parvorder',
    'infraorder',
    'suborder',
    'order',
    'superorder',
    'subcohort',
    'cohort',
    'infraclass',
    'subclass',
    'class',
    'superclass',
    'subphylum',
    'phylum',
    'superphylum',
    'subkingdom',
    'kingdom',
    'superkingdom',
    'root',
]

UNORDERED_RANKS = [
    'clade',
    'isolate',
    'forma_specialis',
    'no_rank',
    'serotype',
    'strain']

RANKS = RANK_ORDER + UNORDERED_RANKS

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
                                 # starts with lower-case char (except root)
                                 r'^(?!root$)[a-z]',
                                 r'^\W+\s+[a-zA-Z]*\d',  # Digit in second word
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
                                 # r'\b[Gg]roup\b',
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


def define_schema(Base):

    class Node(Base):
        __tablename__ = 'nodes'
        tax_id = Column(String, primary_key=True, nullable=False)

        # TODO: temporarily remove foreign key constratint on parent
        # TODO: (creates order depencence during insertion); may need to
        # TODO: add constraint after table has been populated.

        # parent_id = Column(String, ForeignKey('nodes.tax_id'))
        parent_id = Column(String, index=True)
        rank = Column(String, ForeignKey('ranks.rank'))
        embl_code = Column(String)
        division_id = Column(String)
        source_id = Column(Integer, ForeignKey('source.id'))
        is_valid = Column(Boolean, default=True)
        names = relationship('Name')
        ranks = relationship('Rank', back_populates='nodes')
        sources = relationship('Source', back_populates='nodes')

    class Name(Base):
        __tablename__ = 'names'
        # id = Column(Integer, primary_key=True)
        tax_id = Column(String, ForeignKey('nodes.tax_id', ondelete='CASCADE'))
        node = relationship('Node', back_populates='names')
        tax_name = Column(String)
        unique_name = Column(String)
        name_class = Column(String)
        source_id = Column(Integer, ForeignKey('source.id'))
        is_primary = Column(Boolean)
        is_classified = Column(Boolean)
        sources = relationship('Source', back_populates='names')
        __table_args__ = (
            PrimaryKeyConstraint('tax_id', 'tax_name', 'name_class'), {},)

    Index('ix_names_tax_id_is_primary', Name.tax_id, Name.is_primary)

    class Merge(Base):
        __tablename__ = 'merged'
        old_tax_id = Column(String, primary_key=True, index=True)
        new_tax_id = Column(String, ForeignKey(
            'nodes.tax_id', ondelete='CASCADE'))

    class Rank(Base):
        __tablename__ = 'ranks'
        rank = Column(String, primary_key=True)
        height = Column(Integer, unique=True, nullable=False)
        nodes = relationship('Node')

    class Source(Base):
        __tablename__ = 'source'
        id = Column(Integer, primary_key=True)
        name = Column(String, unique=True)
        description = Column(String)
        nodes = relationship('Node')
        names = relationship('Name')


def db_connect(engine, schema=None, clobber=False):
    """Create a connection object to a database. Attempt to establish a
    schema. If there are existing tables, delete them if clobber is
    True and return otherwise. Returns a sqlalchemy engine object.

    """

    if schema is None:
        base = declarative_base()
    else:
        try:
            engine.execute(sqlalchemy.schema.CreateSchema(schema))
        except sqlalchemy.exc.ProgrammingError as err:
            logging.warn(err)
        base = declarative_base(metadata=MetaData(schema=schema))

    define_schema(base)

    if clobber:
        logging.info('Clobbering database tables')
        base.metadata.drop_all(bind=engine)

    logging.info('Creating database tables')
    base.metadata.create_all(bind=engine)

    return base


def read_merged(rows):

    yield ('old_tax_id', 'new_tax_id')
    for row in rows:
        yield tuple(row)


def read_nodes(rows, source_id=1):
    """
    Return an iterator of rows ready to insert into table "nodes".

    * rows - iterator of lists (eg, output from read_archive or read_dmp)
    """

    ncbi_keys = ['tax_id', 'parent_id', 'rank', 'embl_code', 'division_id']
    extra_keys = ['source_id', 'is_valid']
    is_valid = True

    ncbi_cols = len(ncbi_keys)

    rank = ncbi_keys.index('rank')
    parent_id = ncbi_keys.index('parent_id')

    # assumes the first row is the root
    row = next(rows)
    row[rank] = 'root'
    # parent must be None for termination of recursive CTE for
    # calculating lineages
    row[parent_id] = None
    rows = itertools.chain([row], rows)

    yield ncbi_keys + extra_keys

    for row in rows:
        # replace whitespace in "rank" with underscore
        row[rank] = '_'.join(row[rank].split())
        # provide default values for source_id and is_valid
        yield row[:ncbi_cols] + [source_id, is_valid]


def read_names(rows, source_id=1):
    """Return an iterator of rows ready to insert into table
    "names". Adds columns "is_primary" (identifying the primary name
    for each tax_id with a vaule of 1) and "is_classified" (always None).

    * rows - iterator of lists (eg, output from read_archive or read_dmp)
    * unclassified_regex - a compiled re matching "unclassified" names

    From the NCBI docs:

    Taxonomy names file (names.dmp):
        tax_id -- the id of node associated with this name
        name_txt -- name itself
        unique name -- the unique variant of this name if name not unique
        name class -- (synonym, common name, ...)

    """

    ncbi_keys = ['tax_id', 'tax_name', 'unique_name', 'name_class']
    extra_keys = ['source_id', 'is_primary', 'is_classified']

    # is_classified applies to species only; we will set this value
    # later
    is_classified = None

    tax_id = ncbi_keys.index('tax_id')
    tax_name = ncbi_keys.index('tax_name')
    unique_name = ncbi_keys.index('unique_name')
    name_class = ncbi_keys.index('name_class')

    yield ncbi_keys + extra_keys

    for tid, grp in itertools.groupby(rows, itemgetter(tax_id)):
        # confirm that each tax_id has exactly one scientific name
        num_primary = 0
        for r in grp:
            is_primary = r[name_class] == 'scientific name'
            # fix primary key uniqueness violation
            if r[unique_name]:
                r[tax_name] = r[unique_name]
            num_primary += is_primary
            yield (r + [source_id, is_primary, is_classified])

        assert num_primary == 1


class NCBILoader(object):
    def __init__(self, engine, schema=None, ranks=RANKS):
        self.engine = engine
        self.schema = schema
        self.tables = {name: self.prepend_schema(name)
                       for name in ['merged', 'names', 'nodes', 'ranks', 'source']}
        self.ranks = ranks
        self.placeholder = {'pysqlite': '?', 'psycopg2': '%s'}[engine.driver]

    def prepend_schema(self, name):
        """Prepend schema name to 'name' when a schema is specified

        """
        return '.'.join([self.schema, name]) if self.schema else name

    def load_table(self, table, rows, colnames=None, limit=None):
        """Load 'rows' into table 'table'. If 'colnames' is not provided, the
        first element of 'rows' must provide column names.

        """

        conn = self.engine.raw_connection()
        cur = conn.cursor()

        colnames = colnames or next(rows)

        cmd = 'INSERT INTO {table} ({colnames}) VALUES ({placeholders})'.format(
            table=self.tables[table],
            colnames=', '.join(colnames),
            placeholders=', '.join([self.placeholder] * len(colnames)))

        cur.executemany(cmd, itertools.islice(rows, limit))
        conn.commit()

    def load_archive(self, archive):
        """Load data from the zip archive of the NCBI taxonomy.

        """

        # source
        self.load_table(
            'source',
            rows=[('ncbi', DATA_URL)],
            colnames=['name', 'description'],
        )

        conn = self.engine.raw_connection()
        cur = conn.cursor()
        cmd = "select id from {source} where name = 'ncbi'".format(**self.tables)
        cur.execute(cmd)
        source_id = cur.fetchone()[0]

        # ranks
        log.info('loading ranks')
        self.load_table(
            'ranks',
            rows=((rank, i) for i, rank in enumerate(RANKS)),
            colnames=['rank', 'height'],
        )

        # nodes
        logging.info('loading nodes')
        nodes_rows = read_nodes(
            read_archive(archive, 'nodes.dmp'), source_id=source_id)
        self.load_table('nodes', rows=nodes_rows)

        # names
        logging.info('loading names')
        names_rows = read_names(
            read_archive(archive, 'names.dmp'), source_id=source_id)
        self.load_table('names', rows=names_rows)

        # merged
        logging.info('loading merged')
        merged_rows = read_merged(read_archive(archive, 'merged.dmp'))
        self.load_table('merged', rows=merged_rows)

    def set_names_is_classified(self, unclassified_regex=UNCLASSIFIED_REGEX):
        conn = self.engine.raw_connection()
        cur = conn.cursor()

        log.info('retrieving species names')

        tempname = random_name(12)
        tablenames = dict(self.tables, temptab=self.prepend_schema(tempname))

        cmd = """
        SELECT tax_id, tax_name
        FROM {names}
        JOIN {nodes} USING(tax_id)
        WHERE is_primary
        AND rank = 'species'
        """.format(**tablenames)

        cur.execute(cmd)
        species_names = cur.fetchall()
        log.info('found {} species names'.format(len(species_names)))

        log.info('checking for unclassified species names')
        classified_taxids = [(tax_id,) for tax_id, tax_name in species_names
                             if not unclassified_regex.search(tax_name)]

        log.info('{count} names are classified ({pct}%)'.format(
            count=len(classified_taxids),
            pct=round((100.0 * len(classified_taxids)) / len(species_names), 1))
        )

        # insert tax_ids into a temporary table
        cmd = 'CREATE TEMPORARY TABLE "{temptab}" (tax_id text)'.format(**tablenames)
        log.info(cmd)
        cur.execute(cmd)

        log.info('inserting tax_ids into temporary table')
        cmd = 'INSERT INTO "{temptab}" VALUES ({placeholder})'.format(
            placeholder=self.placeholder, **tablenames)
        log.info(cmd)
        cur.executemany(cmd, classified_taxids)

        log.info('creating an index on the temporary table')
        cmd = 'CREATE INDEX ix_{tempname}_tax_id on "{temptab}"(tax_id)'.format(
            tempname=tempname, **tablenames)
        log.info(cmd)
        cur.execute(cmd)

        log.info('updating names.is_classified')

        cmd = """
        UPDATE {names} SET
        is_classified = {placeholder}
        WHERE is_primary
        AND tax_id IN (SELECT tax_id FROM "{temptab}")
        """.format(placeholder=self.placeholder, **tablenames)

        log.info(cmd)
        cur.execute(cmd, (True,))
        conn.commit()

    def set_nodes_is_valid(self):
        conn = self.engine.raw_connection()
        cur = conn.cursor()

        cmd = """
        WITH RECURSIVE descendants AS (
         SELECT tax_id, parent_id, rank
         FROM {nodes}
         WHERE rank = 'species'
         AND tax_id not in (SELECT tax_id FROM {names} WHERE is_classified)
         UNION
         SELECT
         n.tax_id,
         n.parent_id,
         n.rank
         FROM {nodes} n
         INNER JOIN descendants d ON d.tax_id = n.parent_id
        )
        UPDATE {nodes} SET is_valid = {placeholder}
        WHERE tax_id in (SELECT tax_id from descendants)
        """.format(placeholder=self.placeholder, **self.tables)

        log.info('marking invalid nodes')
        log.info(cmd)

        cur.execute(cmd, (False, ))
        conn.commit()


def fetch_data(dest_dir='.', clobber=False, url=DATA_URL):
    """
    Download data from NCBI required to generate local taxonomy
    database. Default url is ncbi.DATA_URL

    * dest_dir - directory in which to save output files (created if necessary).
    * clobber - don't download if False and target of url exists in dest_dir
    * url - url to archive; default is ncbi.DATA_URL

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
        logging.info(fout + ' exists; not downloading')
    else:
        downloaded = True
        logging.info('downloading {} to {}'.format(url, fout))
        request.urlretrieve(url, fout)

    return (fout, downloaded)


def read_archive(archive, fname):
    """Return an iterator of unique rows from a zip archive.

    * archive - path to the zip archive.
    * fname - name of the compressed file within the archive.

    """

    # Note that deduplication here is equivalent to an upsert/ignore,
    # but avoids requirement for a database-specific implementation.

    zfile = zipfile.ZipFile(archive)
    contents = zfile.open(fname, 'r')
    fobj = io.TextIOWrapper(contents)

    seen = set()
    for line in fobj:
        line = line.rstrip('\t|\n')
        if line not in seen:
            yield line.split('\t|\t')
            seen.add(line)

# def read_dmp(fname):
#     seen = set()
#     for line in open(fname, 'rU'):
#         line = line.rstrip('\t|\n')
#         if line not in seen:
#             yield line.split('\t|\t')
#             seen.add(line)
