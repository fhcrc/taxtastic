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
The Taxonomy class defines an object providing an interface to
the taxonomy database.
"""
import csv
import logging
import random
import string

from jinja2 import Template

import sqlalchemy
from sqlalchemy import MetaData, and_, or_
from sqlalchemy.sql import select

log = logging.getLogger(__name__)


class TaxonIntegrityError(Exception):
    '''
    Raised when something in the Taxonomy is not structured correctly
    '''
    pass


class Taxonomy(object):

    def __init__(self, engine, NO_RANK='no_rank',
                 undef_prefix='below_', schema=None):
        """
        The Taxonomy class defines an object providing an interface to
        the taxonomy database.

        * engine - sqlalchemy engine instance providing a connection to a
          database defining the taxonomy
        * NO_RANK - label identifying a taxon without
          a specific rank in the taxonomy.
        * undef_prefix - string prepended to name of parent
          rank to create new labels for undefined ranks.
        * schema - database schema, usually required when using a Postgres db

        Example:
        >>> from sqlalchemy import create_engine
        >>> from taxtastic.taxonomy import Taxonomy
        >>> engine = create_engine(url)
        >>> tax = Taxonomy(engine)

        see http://www.sqlalchemy.org/docs/reference/sqlalchemy/inspector.html
        http://www.sqlalchemy.org/docs/metadata.html#metadata-reflection
        """

        log.debug('using database ' + str(engine.url))

        self.engine = engine
        self.meta = MetaData(schema=schema)
        self.meta.bind = self.engine
        self.meta.reflect()

        schema_prefix = schema + '.' if schema else ''

        self.nodes = self.meta.tables[schema_prefix + 'nodes']
        self.names = self.meta.tables[schema_prefix + 'names']
        self.source = self.meta.tables[schema_prefix + 'source']
        self.merged = self.meta.tables[schema_prefix + 'merged']
        ranks_table = self.meta.tables[schema_prefix + 'ranks']
        ranks = select([ranks_table.c.rank, ranks_table.c.height]).execute().fetchall()
        ranks = sorted(ranks, key=lambda x: int(x[1]))  # sort by height
        self.ranks = [r[0] for r in ranks]  # just the ordered ranks

        if 'taxonomy' in self.meta.tables:
            self.taxonomy = self.meta.tables[schema_prefix + 'taxonomy']
        else:
            self.taxonomy = None

        self.NO_RANK = NO_RANK
        self.undef_prefix = undef_prefix

        # self.cached = {}

        # TODO: can probably remove this check at some point;
        # historically the root node had tax_id == parent_id ==
        # '1', but with a recursive CTE, parent_id must be None for
        # the recursive expression to terminate.
        with self.engine.connect() as con:
            result = con.execute("select parent_id from nodes where rank = 'root'")
            if result.fetchone()[0] is not None:
                raise TaxonIntegrityError('the root node must have parent_id = None')

    def _add_rank(self, rank, parent_rank):
        """
        inserts rank into self.ranks.
        """
        if rank not in self.rankset:
            self.ranks.insert(self.ranks.index(parent_rank) + 1, rank)
        self.rankset = set(self.ranks)

    def _node(self, tax_id):
        """
        Returns parent_id, rank

        FIXME: expand return rank to include custom 'below' ranks built when
               get_lineage is caled
        """
        s = select([self.nodes.c.parent_id, self.nodes.c.rank],
                   self.nodes.c.tax_id == tax_id)
        res = s.execute()
        output = res.fetchone()
        if not output:
            msg = 'value "{}" not found in nodes.tax_id'.format(tax_id)
            raise ValueError(msg)
        else:
            return output  # parent_id, rank

    def primary_from_id(self, tax_id):
        """
        Returns primary taxonomic name associated with tax_id
        """
        s = select([self.names.c.tax_name],
                   and_(self.names.c.tax_id == tax_id,
                        self.names.c.is_primary))
        res = s.execute()
        output = res.fetchone()

        if not output:
            msg = 'value "{}" not found in names.tax_id'.format(tax_id)
            raise ValueError(msg)
        else:
            return output[0]

    def primary_from_name(self, tax_name):
        """
        Return tax_id and primary tax_name corresponding to tax_name.
        """
        names = self.names

        s1 = select([names.c.tax_id, names.c.is_primary],
                    names.c.tax_name == tax_name)

        log.debug(str(s1))

        res = s1.execute().fetchone()
        if res:
            tax_id, is_primary = res
        else:
            msg = '"{}" not found in names.tax_names'.format(tax_name)
            raise ValueError(msg)

        if not is_primary:
            s2 = select([names.c.tax_name],
                        and_(names.c.tax_id == tax_id,
                             names.c.is_primary))
            tax_name = s2.execute().fetchone()[0]

        return tax_id, tax_name, bool(is_primary)

    def _get_merged(self, tax_id):
        """Returns tax_id into which `tax_id` has been merged or `tax_id` of
        not obsolete.

        """

        cmd = """
        SELECT COALESCE(
        (SELECT new_tax_id FROM merged
         WHERE old_tax_id = {x}), {x})
        """.format(x='%s' if self.engine.name == 'postgresql' else '?')

        with self.engine.connect() as con:
            result = con.execute(cmd, (tax_id, tax_id))
            return result.fetchone()[0]

    def _get_lineage(self, tax_id, merge_obsolete=True):
        """Return a list of [(rank, tax_id)] describing the lineage of
        tax_id. If ``merge_obsolete`` is True and ``tax_id`` has been
        replaced, use the corresponding value in table merged.

        """

        try:
            # Be sure we aren't working with an obsolete tax_id
            if merge_obsolete:
                tax_id = self._get_merged(tax_id)

            # Note: joining with ranks seems like a no-op, but for some
            # reason it results in a faster query using sqlite, as well as
            # an ordering from leaf --> root. Might be a better idea to
            # sort explicitly if this is the expected behavior, but it
            # seems like for the most part, the lineage is converted to a
            # dict and the order is irrelevant.
            cmd = """
            WITH RECURSIVE a AS (
             SELECT tax_id, parent_id, rank
              FROM nodes
              WHERE tax_id = {}
            UNION ALL
             SELECT p.tax_id, p.parent_id, p.rank
              FROM a JOIN nodes p ON a.parent_id = p.tax_id
            )
            SELECT a.rank, a.tax_id FROM a
            JOIN ranks using(rank)
            """.format('%s' if self.engine.name == 'postgresql' else '?')

            # reorder so that root is first
            with self.engine.connect() as con:
                result = con.execute(cmd, (tax_id,))
                return result.fetchall()[::-1]
        except sqlalchemy.exc.ResourceClosedError:
            raise ValueError('tax id "{}" not found'.format(tax_id))

    def _get_lineage_table(self, tax_ids, merge_obsolete=True):
        """Return a list of [(rank, tax_id, tax_name)] describing the lineage
        of tax_id. If ``merge_obsolete`` is True and ``tax_id`` has
        been replaced, use the corresponding value in table merged.

        """

        try:
            with self.engine.connect() as con:
                # insert tax_ids into a temporary table
                temptab = ''.join([random.choice(string.ascii_letters)
                                   for n in xrange(12)])
                cmd = 'CREATE TEMPORARY TABLE "{}" (old_tax_id text)'.format(temptab)
                con.execute(cmd)

                # TODO: couldn't find an equivalent of "executemany" - does one exist?
                cmd = 'INSERT INTO "{}" VALUES (?)'.format(temptab)
                for tax_id in tax_ids:
                    con.execute(cmd, tax_id)

                cmd = Template("""
                WITH RECURSIVE a AS (
                 SELECT tax_id as tid, 1 AS ord, tax_id, parent_id, rank
                  FROM nodes
                  WHERE tax_id in (
                  {% if merge_obsolete %}
                  SELECT COALESCE(m.new_tax_id, "{{ temptab }}".old_tax_id)
                    FROM "{{ temptab }}" LEFT JOIN merged m USING(old_tax_id)
                  {% else %}
                  SELECT * from "{{ temptab }}"
                  {% endif %}
                  )
                UNION ALL
                 SELECT a.tid, a.ord + 1, p.tax_id, p.parent_id, p.rank
                  FROM a JOIN nodes p ON a.parent_id = p.tax_id
                )
                SELECT a.tid, a.ord, a.tax_id, a.parent_id, a.rank, tax_name FROM a
                JOIN names using(tax_id)
                WHERE names.is_primary
                ORDER BY tid, ord desc
                """).render(
                    temptab=temptab,
                    merge_obsolete=merge_obsolete,
                )

                result = con.execute(cmd)
                rows = result.fetchall()

                con.execute('DROP TABLE "{}"'.format(temptab))
                return rows

        except sqlalchemy.exc.ResourceClosedError:
            raise ValueError('tax id "{}" not found'.format(tax_id))

    def is_below(self, lower, upper):
        return lower in self.ranks_below(upper)

    def ranks_below(self, rank, depth=None):
        below = []
        try:
            below = self.ranks[self.ranks.index(rank):depth]
        except ValueError as err:
            log.error(err)
        return below

    def synonyms(self, tax_id=None, tax_name=None):
        if not bool(tax_id) ^ bool(tax_name):
            raise ValueError(
                'Exactly one of tax_id and tax_name may be provided.')

        names = self.names

        if tax_name:
            s1 = select([names.c.tax_id], names.c.tax_name == tax_name)
            res = s1.execute().fetchone()

            if res:
                tax_id = res[0]
            else:
                msg = '"{}" not found in names.tax_names'.format(tax_name)
                raise ValueError(msg)

        s = select([names.c.tax_name, names.c.is_primary],
                   names.c.tax_id == tax_id)
        output = s.execute().fetchall()

        if not output:
            raise ValueError('"{}" not found in names.tax_id'.format(tax_id))

        return output

    def lineage(self, tax_id=None, tax_name=None):
        """
        Public method for returning a lineage; includes tax_name and rank
        """

        if not bool(tax_id) ^ bool(tax_name):
            msg = 'Exactly one of tax_id and tax_name may be provided.'
            raise ValueError(msg)

        if tax_name:
            tax_id, primary_name, is_primary = self.primary_from_name(tax_name)
        else:
            primary_name = None

        # assumes stable ordering of lineage from root --> leaf
        lintups = self._get_lineage(tax_id)
        ldict = dict(lintups)

        ldict['tax_id'] = tax_id
        try:
            # parent is second to last element, except for root
            __, ldict['parent_id'] = lintups[-2]
        except IndexError:
            ldict['parent_id'] = None

        ldict['rank'], __ = lintups[-1]  # this taxon is last element in lineage
        ldict['tax_name'] = primary_name or self.primary_from_id(tax_id)

        return ldict

    def add_source(self, name, description=None):
        """
        Attempts to add a row to table "source". Returns (source_id,
        True) if the insert succeeded, (source_id, False) otherwise.
        """

        try:
            result = self.source.insert().execute(name=name,
                                                  description=description)
            source_id, success = result.inserted_primary_key[0], True
        except sqlalchemy.exc.IntegrityError:
            s = select([self.source.c.id], self.source.c.name == name)
            source_id, success = s.execute().fetchone()[0], False

        return source_id, success

    def verify_rank_integrity(self, tax_id, rank, parent_id, children):
        '''
        confirm that for each node the parent ranks and children ranks are coherent
        '''
        def _lower(n1, n2):
            return self.ranks.index(n1) < self.ranks.index(n2)

        if not _lower(rank, self.rank(parent_id)):
            msg = 'New node {} has same or higher rank than parent node {}'
            msg = msg.format(tax_id, 'parent_id')
            raise TaxonIntegrityError(msg)

        for child in children:
            if not _lower(self.rank(child), rank):
                msg = 'Child node {} has same or lower rank as new node {}'
                msg = msg.format(tax_id, child)
                raise TaxonIntegrityError(msg)
        return True

    def add_node(self, tax_id, parent_id, rank, tax_name,
                 children=[], source_id=None, source_name=None):
        """
        Add a node to the taxonomy.
        """

        if rank not in self.ranks:
            msg = 'adding new ranks to taxonomy is not yet supported'
            raise TaxonIntegrityError(msg)

        if not (source_id or source_name):
            raise ValueError(
                'Taxonomy.add_node requires source_id or source_name')

        self.verify_rank_integrity(tax_id, rank, parent_id, children)

        if not source_id:
            source_id, source_is_new = self.add_source(name=source_name)

        self.nodes.insert().execute(tax_id=tax_id,
                                    parent_id=parent_id,
                                    rank=rank,
                                    source_id=source_id)

        self.names.insert().execute(tax_id=tax_id,
                                    tax_name=tax_name,
                                    is_primary=True)

        for child in children:
            ret = self.nodes.update(
                whereclause=self.nodes.c.tax_id == child,
                values={'parent_id': tax_id})
            ret.execute()

        lineage = self.lineage(tax_id)

        log.debug(lineage)

        if self.taxonomy is not None:
            self.taxonomy.insert().execute(
                tax_id=tax_id, parent_id=parent_id, rank=rank, **lineage)

            for child in children:
                c_rank = self.get_rank(child)
                ret = self.taxonomy.update(
                    whereclause=getattr(self.taxonomy.c, c_rank) == child,
                    values={c_rank: tax_id})
                ret.execute()

        return lineage

    def update_node(self, tax_id, **values):
        if all(k not in values for k in ['source_id', 'source_name']):
            msg = 'Taxonomy.update_node requires source_id or source_name: '
            raise ValueError(msg + str(values))
        if 'source_id' not in values:
            source_id, _ = self.add_source(name=values['source_name'])
            values['source_id'] = source_id

        # drop columns not in nodes table
        values = dict(c for c in values.items() if c[0] in self.nodes.c)

        self.nodes.update(
            whereclause=self.nodes.c.tax_id == tax_id,
            values=values).execute()
        lineage = self.lineage(tax_id)
        log.debug(lineage)
        return lineage

    def sibling_of(self, tax_id):
        """Return None or a tax_id of a sibling of *tax_id*.

        If *tax_id* is None, then always returns None. Otherwise,
        returns None if there is no sibling.
        """
        if tax_id is None:
            return None
        parent_id, rank = self._node(tax_id)
        s = select([self.nodes.c.tax_id],
                   and_(self.nodes.c.parent_id == parent_id,
                        self.nodes.c.tax_id != tax_id,
                        self.nodes.c.rank == rank))
        res = s.execute()
        output = res.fetchone()
        if not output:
            msg = 'No sibling of tax_id {} with rank {} found in taxonomy'
            msg = msg.format(tax_id, rank)
            log.warning(msg)
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

    def tax_ids(self):
        '''
        Return all tax_ids in node table
        '''
        fetch = select([self.nodes.c.tax_id]).execute().fetchall()
        ids = [t[0] for t in fetch]
        return ids

    def child_of(self, tax_id):
        """Return None or a tax id of a child of *tax_id*.

        If *tax_id* is None, then always returns None. Otherwise
        returns a child if one exists, else None. The child must have
        a proper rank below that of tax_id (i.e., genus, species, but
        not no_rank or below_below_kingdom).
        """
        if tax_id is None:
            return None
        parent_id, rank = self._node(tax_id)
        s = select([self.nodes.c.tax_id],
                   and_(self.nodes.c.parent_id == tax_id,
                        or_(*[self.nodes.c.rank == r
                              for r in self.ranks_below(rank)])))
        res = s.execute()
        output = res.fetchone()
        if not output:
            msg = ('No children of tax_id {} with '
                   'rank below {} found in database')
            msg = msg.format(tax_id, rank)
            log.warning(msg)
            return None
        else:
            r = output[0]
            assert self.is_ancestor_of(r, tax_id)
            return r

    def children_of(self, tax_id, n):
        if tax_id is None:
            return None
        parent_id, rank = self._node(tax_id)
        s = select([self.nodes.c.tax_id],
                   and_(self.nodes.c.parent_id == tax_id,
                        or_(*[self.nodes.c.rank == r
                              for r in self.ranks_below(rank)]))).limit(n)
        res = s.execute()
        output = res.fetchall()
        if not output:
            return []
        else:
            r = [x[0] for x in output]
            for x in r:
                assert self.is_ancestor_of(x, tax_id)
            return r

    def parent_id(self, tax_id, rank=None):
        parent_id, tax_rank = self._node(tax_id)

        if rank:
            if tax_rank == rank:
                msg = 'tax_id {} already at rank {}, returning None'
                msg = msg.format(tax_id, rank)
                log.warn(msg)

            parent_id = dict(self._get_lineage(parent_id)).get(rank, None)

        return parent_id

    def nary_subtree(self, tax_id, n=2):
        """Return a list of species tax_ids under *tax_id* such that
        node under *tax_id* and above the species has two children.
        """
        if tax_id is None:
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
        if tax_id is None:
            return None
        try:
            parent_id, rank = self._node(tax_id)
        except ValueError:
            return None
        if rank == 'species':
            return tax_id
        else:
            c = self.child_of(tax_id)
            newc = self.species_below(c)
            assert self.is_ancestor_of(newc, tax_id)
            return newc
