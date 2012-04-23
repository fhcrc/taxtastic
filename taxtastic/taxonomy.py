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
import logging
import csv
import itertools

log = logging

import sqlalchemy
from sqlalchemy import MetaData, and_, or_
from sqlalchemy.sql import select

from . import ncbi

class Taxonomy(object):

    def __init__(self, engine, ranks=ncbi.ranks, undefined_rank='no_rank', undef_prefix='below'):
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
        >>> from sqlalchemy import create_engine
        >>> from taxtastic.taxonomy import Taxonomy
        >>> engine = create_engine('sqlite:///%s' % dbname, echo=False)
        >>> tax = Taxonomy(engine)

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
        if tax_id == None:
            return None

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
        """Returns tax_id into which `old_tax_id` has been merged.

        If *old_tax_id* is not obsolete, returns it directly.

        CREATE TABLE merged(
        old_tax_id    TEXT,
        new_tax_id    TEXT REFERENCES nodes(tax_id)
        );
        """

        s = select([self.merged.c.new_tax_id], self.merged.c.old_tax_id == old_tax_id)
        res = s.execute()
        output = res.fetchall() or None

        if output is not None:
            if len(output) > 1:
                raise ValueError('There is more than one value for merged.old_tax_id = "%s"' % old_tax_id)
            else:
                output = output[0][0]
        else:
            output = old_tax_id

        return output

    def _get_lineage(self, tax_id, _level=0, merge_obsolete=True):
        """
        Returns cached lineage from self.cached or recursively builds
        lineage of tax_id until the root node is reached.
        """
        # Be sure we aren't working with an obsolete tax_id
        if merge_obsolete:
            tax_id = self._get_merged(tax_id)

        # Note: indent is referenced through locals() below
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

        new_tax_id = self._get_merged(tax_id)
        if new_tax_id:
            tax_id = new_tax_id

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
            lin = self.cached.values()
        else:
            lin = [self._get_lineage(tax_id) for tax_id in taxa]

        # which ranks are actually represented?
        if full:
            ranks = self.ranks
        else:
            represented = set(itertools.chain.from_iterable(
                    [[node[0] for node in lineage] for lineage in lin]))
            ranks = [r for r in self.ranks if r in represented]

        lineages = [self.lineage(tax_id) for tax_id in taxa]

        fields = ['tax_id', 'parent_id', 'rank', 'tax_name'] + ranks
        writer = csv.DictWriter(csvfile, fieldnames=fields,
                                extrasaction='ignore', quoting=csv.QUOTE_NONNUMERIC)

        # header row
        writer.writeheader()

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

    def add_node(self, tax_id, parent_id, rank, tax_name, children = None, source_id=None, source_name=None, **kwargs):

        """
        Add a node to the taxonomy.
        """

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

        if children:
            for child in children:
                ret = self.nodes.update(self.nodes.c.tax_id == child, {'parent_id':tax_id})
                ret.execute()

        lineage = self.lineage(tax_id)

        log.debug(lineage)
        return lineage

    def sibling_of(self, tax_id):
        """Return None or a tax_id of a sibling of *tax_id*.

        If *tax_id* is None, then always returns None. Otherwise,
        returns None if there is no sibling.
        """
        if tax_id == None:
            return None
        parent_id, rank = self._node(tax_id)
        s = select([self.nodes.c.tax_id],
                   and_(self.nodes.c.parent_id == parent_id,
                        self.nodes.c.tax_id != tax_id,
                        self.nodes.c.rank == rank))
        res = s.execute()
        output = res.fetchone()
        if not output:
            log.warning('No sibling of tax_id %s with rank %s found in taxonomy' % (tax_id, rank))
            return None
        else:
            return output[0]

    def is_ancestor_of(self, node, ancestor):
        if node is None or ancestor is None:
            return False
        l = self.lineage(node)
        return ancestor in l.values()

    def rank(self, tax_id):
        if tax_id is None:
            return None
        else:
            q = self._node(tax_id)
            if q:
                return self._node(tax_id)[1]
            else:
                return None

    def child_of(self, tax_id):
        """Return None or a tax id of a child of *tax_id*.

        If *tax_id* is None, then always returns None. Otherwise
        returns a child if one exists, else None. The child must have
        a proper rank below that of tax_id (i.e., genus, species, but
        not no_rank or below_below_kingdom).
        """
        if tax_id == None:
            return None
        parent_id, rank = self._node(tax_id)
        s = select([self.nodes.c.tax_id],
                   and_(self.nodes.c.parent_id == tax_id,
                        or_(*[self.nodes.c.rank == r for r in ranks_below(rank)])))
        res = s.execute()
        output = res.fetchone()
        if not output:
            log.warning("No children of tax_id %s with rank below %s found in database" % (tax_id, rank))
            return None
        else:
            r = output[0]
            assert self.is_ancestor_of(r, tax_id)
            return r

    def children_of(self, tax_id, n):
        if tax_id == None:
            return None
        parent_id, rank = self._node(tax_id)
        s = select([self.nodes.c.tax_id],
                   and_(self.nodes.c.parent_id == tax_id,
                        or_(*[self.nodes.c.rank == r for r in ranks_below(rank)]))).limit(n)
        res = s.execute()
        output = res.fetchall()
        if not output:
            return []
        else:
            r = [x[0] for x in output]
            for x in r:
                assert self.is_ancestor_of(x, tax_id)
            return r

    def parent_id(self, tax_id):
        if tax_id is None:
            return None
        else:
            q = self._node(tax_id)
            if q:
                return q[0]
            else:
                return None


    def nary_subtree(self, tax_id, n=2):
        """Return a list of species tax_ids under *tax_id* such that
        node under *tax_id* and above the species has two children.
        """
        if tax_id == None:
            return None
        parent_id, rank = self._node(tax_id)
        if rank == 'species':
            return [tax_id]
        else:
            children = self.children_of(tax_id, 2)
            species_taxids = []
            for t in children:
                species_taxids.extend(self.nary_subtree(t, n))
            return species_taxids

    def species_below(self, tax_id):
        if tax_id == None:
            return None
        try:
            parent_id, rank = self._node(tax_id)
        except KeyError:
            return None
        if rank == 'species':
            return tax_id
        else:
            c = self.child_of(tax_id)
            newc = self.species_below(c)
            assert self.is_ancestor_of(newc, tax_id)
            return newc

ranks = ['species', 'genus', 'family', 'order', 'class', 'phylum', 'kingdom']

def is_below(lower, upper):
    try:
        lindex = ranks.index(lower)
        uindex = ranks.index(upper)
        return lindex < uindex
    except:
        return False

def ranks_below(rank):
    try:
        idx = ranks.index(rank)
        return ranks[:idx]
    except:
        return []

def rank_below(rank):
    return {'kingdom': 'phylum',
            'phylum': 'class',
            'class': 'order',
            'order': 'family',
            'family': 'genus',
            'genus': 'species'}[rank]

