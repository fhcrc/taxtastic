import logging
import csv
import itertools
import pprint

log = logging

import sqlalchemy
from sqlalchemy import MetaData, create_engine, and_
from sqlalchemy.sql import select

class Taxonomy(object):

    def __init__(self, engine, ranks, undefined_rank='no_rank', undef_prefix='below'):
        """
        The Taxonomy class defines an object providing an interface to
        the taxonomy database.

        * engine - sqlalchemy engine instance providing a connection to a
          database defining the taxonomy
        * ranks - list of rank names, root first
        * undefined_rank - label identifying a taxon without
          a specific rank in the taxonomy.
        * undef_prefix - string prepended to name of parent
          rank to create new labels for undefined ranks.

        Example:
        > from sqlalchemy import create_engine
        > from taxonomy import Taxonomy, ncbi
        > engine = create_engine('sqlite:///%s' % dbname, echo=False)
        > tax = Taxonomy(engine, taxonomy.ncbi.ranks)

          """

        # TODO: should ranks be defined in a table in the database?
        # TODO: assertions to check for database components

        # see http://www.sqlalchemy.org/docs/reference/sqlalchemy/inspector.html
        # http://www.sqlalchemy.org/docs/metadata.html#metadata-reflection

        log.debug('using database %s' % engine.url)

        self.engine = engine
        self.meta = MetaData()
        self.meta.bind = self.engine
        self.meta.reflect()

        self.nodes = self.meta.tables['nodes']
        self.names = self.meta.tables['names']
        self.source = self.meta.tables['source']
        self.merged = self.meta.tables['merged']

        self.ranks = ranks
        self.rankset = set(self.ranks)

        # keys: tax_id
        # vals: lineage represented as a list of tuples: (rank, tax_id)
        self.cached = {}

        # keys: tax_id
        # vals: lineage represented as a dict of {rank:tax_id}
        # self.taxa = {}

        self.undefined_rank = undefined_rank
        self.undef_prefix = undef_prefix

    def _add_rank(self, rank, parent_rank):
        """
        inserts rank into self.ranks.
        """

        if rank not in self.rankset:
            self.ranks.insert(self.ranks.index(parent_rank) + 1, rank)
        self.rankset = set(self.ranks)

    def _node(self, tax_id):
        """
        Returns parent, rank
        """

        s = select([self.nodes.c.parent_id, self.nodes.c.rank],
                   self.nodes.c.tax_id == tax_id)
        res = s.execute()
        output = res.fetchone()

        if not output:
            raise KeyError('value "%s" not found in nodes.tax_id' % tax_id)

        # parent_id, rank
        return output

    def primary_from_id(self, tax_id):
        """
        Returns primary taxonomic name associated with tax_id
        """

        s = select([self.names.c.tax_name],
                   and_(self.names.c.tax_id == tax_id, self.names.c.is_primary == 1))
        res = s.execute()
        output = res.fetchone()

        if not output:
            raise KeyError('value "%s" not found in names.tax_id' % tax_id)
        else:
            return output[0]

    def primary_from_name(self, tax_name):
        """
        Return tax_id and primary tax_name corresponding to tax_name.
        """

        names = self.names

        s1 = select([names.c.tax_id, names.c.is_primary], names.c.tax_name == tax_name)

        res = s1.execute().fetchone()
        if res:
            tax_id, is_primary = res
        else:
            raise KeyError('"%s" not found in names.tax_names' % tax_name)

        if not is_primary:
            s2 = select([names.c.tax_name],
                        and_(names.c.tax_id == tax_id, names.c.is_primary == 1))
            tax_name = s2.execute().fetchone()[0]

        return tax_id, tax_name, bool(is_primary)

    def _get_merged(self, old_tax_id):
        """
        Returns tax_id into which `old_tax_id` has been merged.

        CREATE TABLE merged(
        old_tax_id    TEXT,
        new_tax_id    TEXT REFERENCES nodes(tax_id)
        );
        """

        merged = self.merged
        s = select([self.merged.c.new_tax_id], self.merged.c.old_tax_id == old_tax_id)
        res = s.execute()
        output = res.fetchall() or None

        if output is not None:
            if len(output) > 1:
                raise ValueError('There is more than one value for merged.old_tax_id = "%s"' % old_tax_id)
            else:
                output = output[0][0]
            
        return output
    
    def _get_lineage(self, tax_id, _level=0):
        """
        Returns cached lineage from self.cached or recursively builds
        lineage of tax_id until the root node is reached.
        """

        indent = '.'*_level

        undefined = self.undefined_rank
        prefix = self.undef_prefix+'_'

        lineage = self.cached.get(tax_id)

        if lineage:
            log.debug('%(indent)s tax_id "%(tax_id)s" is cached' % locals())
        else:
            log.debug('%(indent)s reconstructing lineage of tax_id "%(tax_id)s"' % locals())
            parent_id, rank = self._node(tax_id)
            lineage = [(rank, tax_id)]

            # recursively add parent_ids until we reach the root
            if parent_id != tax_id:
                lineage = self._get_lineage(parent_id, _level+1) + lineage

            # now that we've reached the root, rename any undefined ranks
            _parent_rank, _parent_id = None, None
            for i, node in enumerate(lineage):
                _rank, _tax_id = node

                if _rank == undefined:
                    _rank = prefix + _parent_rank
                    self._add_rank(_rank, _parent_rank)

                    lineage[i] = (_rank, _tax_id)
                    self.cached[_tax_id] = lineage
                    log.debug('renamed undefined rank to %(_rank)s in element %(i)s of lineage of %(tax_id)s' \
                                  % locals())

                _parent_rank, _parent_id = _rank, _tax_id

            self.cached[tax_id] = lineage

        return lineage

    def synonyms(self, tax_id=None, tax_name=None):
        if not bool(tax_id) ^ bool(tax_name):
            raise ValueError('Exactly one of tax_id and tax_name may be provided.')

        names = self.names

        if tax_name:
            s1 = select([names.c.tax_id], names.c.tax_name == tax_name)
            res = s1.execute().fetchone()

            if res:
                tax_id = res[0]
            else:
                raise KeyError('"%s" not found in names.tax_names' % tax_name)

        s = select([names.c.tax_name, names.c.is_primary],
                   names.c.tax_id == tax_id)
        output = s.execute().fetchall()

        if not output:
            raise KeyError('"%s" not found in names.tax_id' % tax_id)

        return output


    def lineage(self, tax_id=None, tax_name=None):
        """
        Public method for returning a lineage; includes tax_name and rank
        """

        if not bool(tax_id) ^ bool(tax_name):
            raise ValueError('Exactly one of tax_id and tax_name may be provided.')

        if tax_name:
            tax_id, primary_name, is_primary = self.primary_from_name(tax_name)

        ldict = dict(self._get_lineage(tax_id))

        ldict['tax_id'] = tax_id
        ldict['parent_id'], _ = self._node(tax_id)
        ldict['rank'] = self.cached[tax_id][-1][0]
        ldict['tax_name'] = self.primary_from_id(tax_id)

        return ldict

    def write_table(self, taxa=None, csvfile=None, full=False):
        """
        Represent the currently defined taxonomic lineages as a rectangular
        array with columns named "tax_id","rank","tax_name", followed
        by a column for each rank proceeding from the root to the more
        specific ranks.

         * taxa - list of taxids to include in the output; if none are
           provided, use self.cached.keys() (ie, those taxa loaded into the cache).
         * csvfile - an open file-like object (see "csvfile" argument to csv.writer)
         * full - if True (the default), includes a column for each rank in self.ranks;
           otherwise, omits ranks (columns) the are undefined for all taxa.
        """

        if not taxa:
            taxa = self.cached.keys()

        # which ranks are actually represented?
        if full:
            ranks = self.ranks
        else:
            represented = set(itertools.chain.from_iterable(
                    [[node[0] for node in lineage] for lineage in self.cached.values()])
            )
            ranks = [r for r in self.ranks if r in represented]

        fields = ['tax_id','parent_id','rank','tax_name'] + ranks
        writer = csv.DictWriter(csvfile, fieldnames=fields,
                                extrasaction='ignore', quoting=csv.QUOTE_NONNUMERIC)

        # header row
        writer.writerow(dict(zip(fields, fields)))
        lineages = [self.lineage(tax_id) for tax_id in taxa]

        for lin in sorted(lineages, key=lambda x: (ranks.index(x['rank']), x['tax_name'])):
             writer.writerow(lin)

    def add_source(self, name, description=None):
        """
        Attempts to add a row to table "source". Returns (source_id,
        True) if the insert succeeded, (source_id, False) otherwise.
        """

        try:
            result = self.source.insert().execute(name = name, description = description)
            source_id, success = result.inserted_primary_key[0], True
        except sqlalchemy.exc.IntegrityError:
            s = select([self.source.c.id], self.source.c.name == name)
            source_id, success = s.execute().fetchone()[0], False

        return source_id, success

    def add_node(self, tax_id, parent_id, rank, tax_name, source_id=None, source_name=None, **kwargs):

        if not (source_id or source_name):
            raise ValueError('Taxonomy.add_node requires source_id or source_name')

        if not source_id:
            source_id, source_is_new = self.add_source(name=source_name)

        result = self.nodes.insert().execute(tax_id = tax_id,
                                             parent_id = parent_id,
                                             rank = rank,
                                             source_id = source_id)

        result = self.names.insert().execute(tax_id = tax_id,
                                             tax_name = tax_name,
                                             is_primary = 1)

        lineage = self.lineage(tax_id)

        log.debug(lineage)
        return lineage

