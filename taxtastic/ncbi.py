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
import pandas
import re
import urllib
import zipfile

import sqlalchemy
from sqlalchemy import (Column, Integer, String, Boolean,
                        ForeignKey, Index, MetaData)
from sqlalchemy.orm import relationship
from sqlalchemy.exc import IntegrityError
from sqlalchemy.ext.declarative import declarative_base

DATA_URL = 'ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdmp.zip'

# For rank order: https://en.wikipedia.org/wiki/Taxonomic_rank
RANKS = [
    'forma',
    'varietas',
    'subspecies',
    'species',
    'species_subgroup',
    'species_group',
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
    'root'
]


# Components of a regex to apply to all names. Names matching this regex are
# marked as invalid.
UNCLASSIFIED_REGEX_COMPONENTS = [r'-like\b',
                                 r'\bcrenarchaeote\b',
                                 r'\bculture\b',
                                 r'\bchimeric\b',
                                 # r'\b[cC]lone',
                                 r'\bmethanogen(ic)?\b',
                                 r'\bplanktonic\b',
                                 r'\bplanctomycete\b',
                                 r'\bsymbiote\b',
                                 r'\btransconjugant\b',
                                 # starts with lower-case char (except root)
                                 r'^(?!root$)[a-z]',
                                 # r'^\W+\s+[a-zA-Z]*\d',  # Digit in second word
                                 r'\d\d',
                                 # r'atypical',
                                 # r'^cf\.',
                                 # r'acidophile',
                                 # r'\bactinobacterium\b',
                                 # r'aerobic',
                                 r'.+\b[Al]g(um|a)\b',
                                 # r'\b[Bb]acteri(um|al)\b',
                                 # r'.+\b[Bb]acteria\b',
                                 r'Barophile',
                                 # r'cyanobacterium',
                                 # r'Chloroplast',
                                 r'Cloning',
                                 # r'\bclone\b',
                                 r'cluster',
                                 r'^diazotroph',
                                 # r'\bcoccus\b',
                                r'archaeon',
                                 r'-containing',
                                 r'epibiont',
                                 # 'et al',
                                 r'environmental samples',
                                 # r'eubacterium',
                                 # r'\b[Gg]roup\b',
                                 # r'halophilic',
                                 # r'hydrothermal\b',
                                 # r'isolate',
                                 # r'\bmarine\b',
                                 r'methanotroph',
                                 # r'microorganism',
                                 r'mollicute',
                                 # r'pathogen',
                                 r'[Pp]hytoplasma',
                                 # r'proteobacterium',
                                 r'putative',
                                 # r'\bsp\.',
                                 # r'species',
                                 # r'spirochete',
                                 # r'str\.',
                                 # r'strain',
                                 # r'symbiont',
                                 # r'\b[Tt]axon\b',
                                 # r'unicellular',
                                 # r'uncultured',
                                 # r'unclassified',
                                 # r'unidentified',
                                 # r'unknown',
                                 r'vector\b',
                                 r'vent\b',
                                 ]

# provides criteria for defining matching tax_ids as "unclassified"
UNCLASSIFIED_REGEX = re.compile('|'.join(UNCLASSIFIED_REGEX_COMPONENTS))


def define_schema(Base):
    class Node(Base):
        __tablename__ = 'nodes'
        tax_id = Column(String, primary_key=True, nullable=False)
        parent_id = Column(String, ForeignKey('nodes.tax_id'))
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
        id = Column(Integer, primary_key=True)
        tax_id = Column(String, ForeignKey('nodes.tax_id', ondelete='CASCADE'))
        node = relationship('Node', back_populates='names')
        tax_name = Column(String)
        unique_name = Column(String)
        name_class = Column(String)
        source_id = Column(Integer, ForeignKey('source.id'))
        is_primary = Column(Boolean)
        is_classified = Column(Boolean)
        sources = relationship('Source', back_populates='names')

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
        no_rank = Column(Boolean)
        nodes = relationship('Node')

    class Source(Base):
        __tablename__ = 'source'
        id = Column(Integer, primary_key=True)
        name = Column(String, unique=True)
        description = Column(String)
        nodes = relationship('Node')
        names = relationship('Name')


def db_connect(engine, schema=None, clobber=False):
    """
    Create a connection object to a database. Attempt to establish a
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


def db_load(engine, archive, schema=None):
    """
    Load data from zip archive into database identified by con. Data
    is not loaded if target tables already contain data.
    """

    # source
    source = pandas.DataFrame(
        data={'name': 'ncbi', 'description': DATA_URL}, index=[1])
    source.index.name = 'id'

    # names
    logging.info("Reading names from archive")
    rows = read_names(
        rows=read_archive(archive, 'names.dmp'),
        unclassified_regex=UNCLASSIFIED_REGEX)
    names = pandas.DataFrame(rows)
    names['source_id'] = 1

    assert_primaries(names)  # this should always exist

    # nodes
    logging.info("Reading nodes from archive")
    rows = read_nodes(
        rows=read_archive(archive, 'nodes.dmp'),
        ncbi_source_id=1)
    nodes = pandas.DataFrame(rows).set_index('tax_id')

    logging.info("Marking nodes validity based on primary name")
    nodes = mark_is_valid(nodes, names)

    # add parent rank column for rank adjustments
    nodes = nodes.join(nodes['rank'], on='parent_id', rsuffix='_parent')

    logging.info('Adjusting taxons with same rank as parent')
    nodes = adjust_same_ranks(nodes)

    logging.info('Expanding `no_rank` taxons')
    nodes, ranks = adjust_node_ranks(nodes, RANKS[:])

    # set node ranks as a sortable categories, forma < ... < root
    nodes['rank'] = nodes['rank'].astype('category', categories=ranks, ordered=True)
    nodes['rank_parent'] = nodes['rank_parent'].astype(
        'category', categories=ranks, ordered=True)

    logging.info('Confirming tax tree rank integrity')
    assert_integrity(nodes, ranks)

    logging.info('Marking species subtree validity')
    nodes = mark_valid_subtrees(nodes)

    # merged
    rows = read_archive(archive, 'merged.dmp')
    rows = (dict(zip(['old_tax_id', 'new_tax_id'], row)) for row in rows)
    merged = pandas.DataFrame(rows)

    # ## prepare nodes
    nodes = nodes.drop('rank_parent', axis=1)
    # source_id = 1
    nodes['source_id'] = 1
    # to avoid parent_id foreign key constrant Integrity errors
    nodes = nodes.sort_values('rank', ascending=False)

    logging.info('Inserting source')
    source.to_sql(
        'source', engine,
        schema=schema,
        if_exists='append')

    '''lowest ranks first (forma, species, etc)
    highest ranks last (root, superkingdom, etc)'''
    logging.info('Inserting ranks')
    ranks_df = pandas.DataFrame(data=ranks, columns=['rank'])
    ranks_df['no_rank'] = ranks_df['rank'].apply(lambda x: 'below' in x)
    ranks_df.index.name = 'height'
    ranks_df.to_sql(
        'ranks', engine,
        schema=schema,
        if_exists='append')

    logging.info("Inserting nodes")
    nodes.to_sql(
        'nodes', engine,
        schema=schema,
        if_exists='append')

    logging.info("Inserting names")
    names.to_sql(
        'names', engine,
        schema=schema,
        if_exists='append',
        index=False)

    logging.info("Inserting merged")
    merged.to_sql(
        'merged', engine,
        schema=schema,
        if_exists='append',
        index=False)


def adjust_node_ranks(df, ranks):
    '''
    replace no_ranks with below_ of parent rank
    '''
    ranks = ranks[::-1]  # start at root, work down
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
    return df, ranks[::-1]  # return ranks to input order


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


def assert_integrity(nodes, ranks):
    bad_nodes = nodes[nodes['rank_parent'] < nodes['rank']]
    if not bad_nodes.empty:
        logging.error(bad_nodes)
        raise IntegrityError('some node ranks above parent ranks')


def assert_primaries(names):
    '''
    Search tax_id groups for missing primary names,
    if exist raise IntegrityException
    '''
    err = False
    for i, df in names.groupby(by='tax_id'):
        if not df['is_primary'].any():
            err = True
            logging.error('tax_id {} missing primary name'.format(i))

    if err:
        raise IntegrityError('taxon groups missing primary name')


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
        urllib.urlretrieve(url, fout)

    return (fout, downloaded)


def get_subtrees(nodes, to_mark):
    '''
    to_mark - Boolean Series of nodes for selection

    Continuously grab children using parent_id and OR results on nodes table
    to identify subtrees
    '''
    mark_children = to_mark
    while mark_children.any():
        mark_children = nodes['parent_id'].isin(nodes[mark_children].index)
        to_mark |= mark_children
    return to_mark


def mark_is_valid(nodes, names):
    """
    After applying a ``regex`` to primary names associated with tax_ids,
    marking those that match as invalid.
    """
    nodes = nodes.drop('is_valid', axis=1)
    nodes = nodes.merge(
        names[names['is_primary']][['is_classified', 'tax_id']],
        right_on='tax_id',
        left_index=True,
        how='left')
    nodes = nodes.rename(columns={'is_classified': 'is_valid'})
    return nodes.set_index('tax_id')


def mark_valid_subtrees(nodes):
    subtree_msg = '{} species subtree nodes is_valid={}'

    valid_subtrees = ((nodes['rank'] == 'species') & nodes['is_valid'])
    valid_subtrees = get_subtrees(nodes, valid_subtrees)
    nodes.loc[valid_subtrees, 'is_valid'] = True
    logging.info(subtree_msg.format(valid_subtrees.sum(), True))

    # mark false
    invalid_subtrees = ((nodes['rank'] == 'species') & ~nodes['is_valid'])
    invalid_subtrees = get_subtrees(nodes, invalid_subtrees)
    nodes.loc[invalid_subtrees, 'is_valid'] = False
    logging.info(subtree_msg.format(invalid_subtrees.sum(), False))

    # unclassfied bacteria, tax_id 2323
    unclassified_bacteria = get_subtrees(nodes, nodes.index == '2323')
    logging.info(subtree_msg.format(unclassified_bacteria.sum(), False))
    nodes.loc[unclassified_bacteria, 'is_valid'] = False

    return nodes


def _isame_ancestor_rank(df):
    '''
    return boolean series whether row has same rank as parent
    '''
    return ((df['rank'] != 'root') &
            (df['rank_parent'] != 'no_rank') &
            (df['rank_parent'] == df['rank']))


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

    keys = ['tax_id', 'tax_name', 'unique_name', 'name_class']
    idx = dict((k, i) for i, k in enumerate(keys))
    tax_name, unique_name, name_class = \
        [idx[k] for k in ['tax_name', 'unique_name', 'name_class']]

    def _is_primary(row):
        """
        Defines a name as "primary," meaning that other names associated
        with the this tax_id should be considered synonyms.
        """

        if row[name_class] != 'scientific name':
            is_primary = False
        elif not row[unique_name]:
            is_primary = True
        elif row[tax_name] == row[unique_name].split('<')[0].strip():
            is_primary = True
        else:
            is_primary = False

        return is_primary

    if unclassified_regex:
        def _is_classified(row):
            """
            Return 1 if tax_name element of `row` matches
            unclassified_regex, 0 otherwise. Search no more than the
            first two whitespace-delimited words.
            """
            tn = row[tax_name]
            return not unclassified_regex.search(tn)
    else:
        def _is_classified(row):
            return None

    # appends additional field is_primary
    colnames = keys + ['is_primary', 'is_classified']
    for r in rows:
        row = dict(zip(colnames, r + [_is_primary(r), _is_classified(r)]))
        yield row


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
        row['is_valid'] = True  # default is_valid
        yield row
