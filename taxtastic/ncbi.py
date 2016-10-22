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
import pandas
import re
import urllib
import zipfile

import sqlalchemy
from sqlalchemy import Column, Integer, String, Boolean, ForeignKey, Index
from sqlalchemy.orm import relationship
from sqlalchemy.ext.declarative import declarative_base

from errors import IntegrityError

log = logging

Base = declarative_base()

DATA_URL = 'ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdmp.zip'


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

    tax_id = Column(String, ForeignKey(
        'nodes.tax_id', ondelete='CASCADE'), index=True)
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
    new_tax_id = Column(String, ForeignKey(
        'nodes.tax_id', ondelete='CASCADE'), index=True)
    merged_node = relationship('Node', backref='merged_ids')


class Source(Base):
    __tablename__ = 'source'

    id = Column(Integer, primary_key=True)
    name = Column(String, unique=True)
    description = Column(String)


RANKS = [
    'root',
    'superkingdom',
    'kingdom',
    'subkingdom',
    'superphylum',
    'phylum',
    'subphylum',
    'superclass',
    'class',
    'subclass',
    'infraclass',
    'superorder',
    'order',
    'suborder',
    'infraorder',
    'parvorder',
    'superfamily',
    'family',
    'subfamily',
    'tribe',
    'subtribe',
    'genus',
    'subgenus',
    'species_group',
    'species_subgroup',
    'species',
    'subspecies',
    'varietas',
    'forma'
]

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
                                 r'^[a-z]',  # starts with lower-case character
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


def db_load(engine, archive):
    """
    Load data from zip archive into database identified by con. Data
    is not loaded if target tables already contain data.
    """

    try:
        # nodes
        rows = read_nodes(
            rows=read_archive(archive, 'nodes.dmp'),
            ncbi_source_id=1)
        nodes = pandas.DataFrame(rows).set_index('tax_id')
        nodes = nodes.join(nodes['rank'], on='parent_id', rsuffix='_parent')

        log.info('Adjusting taxons with same rank as parent')
        nodes = adjust_same_ranks(nodes)

        log.info('Expanding `no_rank` taxons')
        nodes, ranks = adjust_node_ranks(nodes, RANKS[:])

        log.info('Confirming tax tree rank integrity')
        nodes['rank'] = nodes['rank'].astype('category', categories=ranks, ordered=True)
        nodes['rank_parent'] = nodes['rank_parent'].astype(
            'category', categories=ranks, ordered=True)
        bad_nodes = nodes[nodes['rank_parent'] < nodes['rank']]
        if not bad_nodes.empty:
            print(bad_nodes)
            raise IntegrityError('some node ranks above parent ranks')

        log.info('Inserting ranks')
        # reverse list so lowest is first and gets the 0 index
        pandas.Series(ranks, dtype=str, name='rank').to_sql(
            'ranks', engine,
            index=False,
            flavor='sqlite',
            if_exists='append')

        log.info("Inserting nodes")
        nodes.drop('rank_parent', axis=1).to_sql(
            'nodes', engine,
            flavor='sqlite',
            if_exists='append')

        # names
        rows = read_names(
            rows=read_archive(archive, 'names.dmp'),
            unclassified_regex=UNCLASSIFIED_REGEX)
        names = pandas.DataFrame(rows)
        log.info("Inserting names")
        names.to_sql(
            'names', engine,
            flavor='sqlite',
            if_exists='append',
            index=False)

        # merged
        rows = read_archive(archive, 'merged.dmp')
        rows = (dict(zip(['old_tax_id', 'new_tax_id'], row)) for row in rows)
        merged = pandas.DataFrame(rows)
        log.info("Inserting merged")
        merged.to_sql(
            'merged', engine,
            flavor='sqlite',
            if_exists='append',
            index=False)

        fix_missing_primary(engine)

        # Mark names as valid/invalid
        log.info("Marking nodes validity based on primary name")
        mark_is_valid(engine)

        log.info("Updating subtree validity")
        update_subtree_validity(engine)

    except sqlalchemy.exc.IntegrityError as err:
        raise IntegrityError(err)


def _isame_ancestor_rank(df):
    '''
    return boolean series whether row has same rank as parent
    '''
    return ((df['rank'] != 'root') &
            (df['rank_parent'] != 'no_rank') &
            (df['rank_parent'] == df['rank']))


def adjust_same_ranks(df):
    '''
    reset bump parent_id to parent of parent for rows where rank is the same
    as the parent rank
    '''
    same_rank = _isame_ancestor_rank(df)
    while _isame_ancestor_rank(df).any():
        parents = df[same_rank].join(df, on='parent_id', rsuffix='_new')
        df.loc[same_rank, 'parent_id'] = parents['parent_id_new']
        df.loc[same_rank, 'rank_parent'] = parents['rank_parent_new']
        same_rank = _isame_ancestor_rank(df)

    return df


def adjust_node_ranks(df, ranks):
    '''
    replace no_ranks with below_ of parent rank
    '''
    for i, rank in enumerate(ranks):
        no_rank = ((df['rank_parent'] == rank) &
                   (df['rank'] == 'no_rank'))
        if no_rank.any():
            below_rank = 'below_' + rank
            iat_rank = df[no_rank].index
            df.at[iat_rank, 'rank'] = below_rank
            # update `rank_parent` columns
            below_children = df['parent_id'].isin(iat_rank)
            df.loc[below_children, 'rank_parent'] = below_rank
            ranks.insert(i + 1, below_rank)

    # remove non-existent ranks
    node_ranks = set(df['rank'].tolist())
    ranks = [r for r in ranks if r in node_ranks]
    ranks = ranks[::-1]  # reverse order so smallest rank is first

    return df, ranks


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
        if tax_ids:
            log.warn("%d records lack primary names", len(tax_ids))

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
            log.warn("No primary name for tax_id %s. Arbitrarily using %s[%s].",
                     tax_id, tax_name, name_class)
            cursor.execute("""UPDATE names
                SET is_primary = 1
                WHERE tax_id = ? AND tax_name = ? AND
                    unique_name = ? AND name_class = ?""",
                           [tax_id, tax_name, unique_name, name_class])


def mark_is_valid(engine):
    """
    Apply ``regex`` to primary names associated with tax_ids, marking those
    that match as invalid.
    """
    with engine.begin() as connection:
        sql = ('UPDATE nodes SET is_valid = '
               '(SELECT is_classified '
               'FROM names '
               'WHERE names.tax_id = nodes.tax_id AND names.is_primary = 1)')
        log.debug(sql)
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
        log.info("Marking %d subtrees as is_valid=%s",
                 len(to_mark), is_valid)
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
        log.info("marking subtrees for unclassified Bacteria invalid")
        for i, in result:
            mark_subtrees(conn, [i], 0)


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


def read_nodes(rows, ncbi_source_id):
    """
    Return an iterator of rows ready to insert into table "nodes".

    * rows - iterator of lists (eg, output from read_archive or read_dmp)
    * root_name - string identifying the root node (replaces NCBI's default).
    """

    keys = 'tax_id parent_id rank embl_code division_id'.split()
    idx = dict((k, i) for i, k in enumerate(keys))
    rank = idx['rank']

    # assume the first row is the root
    row = rows.next()
    row[rank] = 'root'
    rows = itertools.chain([row], rows)

    colnames = keys + ['source_id'] + ['is_valid']
    for r in rows:
        row = dict(zip(colnames, r))
        assert len(row) == len(colnames)

        # replace whitespace in "rank" with underscore
        row['rank'] = '_'.join(row['rank'].split())
        row['is_valid'] = 1  # default is_valid
        yield row


def read_names(rows, unclassified_regex=None):
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
    idx = dict((k, i) for i, k in enumerate(keys))
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
        def _is_classified(row):
            return None

    # appends additional field is_primary
    colnames = keys + ['is_primary', 'is_classified']
    for r in rows:
        row = dict(zip(colnames, r + [_is_primary(r), _is_classified(r)]))
        yield row
