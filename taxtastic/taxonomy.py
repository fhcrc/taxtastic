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

import logging

from jinja2 import Template

import sqlalchemy
from sqlalchemy import MetaData, and_, or_
from sqlalchemy.sql import select
from sqlalchemy.exc import IntegrityError
from sqlalchemy.orm.session import sessionmaker

from taxtastic.ncbi import UNORDERED_RANKS
from taxtastic.utils import random_name

log = logging.getLogger(__name__)


class TaxonIntegrityError(Exception):
    '''
    Raised when something in the Taxonomy is not structured correctly
    '''
    pass


class Taxonomy(object):

    def __init__(self, engine, schema=None):

        """
        The Taxonomy class defines an object providing an interface to
        the taxonomy database.

        * engine - sqlalchemy engine instance providing a connection to a
          database defining the taxonomy
        * NO_RANK - label identifying a taxon without
          a specific rank in the taxonomy.
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

        self.schema = schema

        # TODO: table names should probably be provided in a dict as a
        # single attribute self.tablenames. This will avoid name
        # collisions (eg with self.ranks) and allow the idiom
        # cmd = 'select * from {<table>}'.format(**self.tablenames)
        self.nodes = self._get_table('nodes')
        self.names = self._get_table('names')
        self.source = self._get_table('source')
        self.merged = self._get_table('merged')
        ranks_table = self._get_table('ranks')
        self.ranks_table = ranks_table

        ranks = select([ranks_table.c.rank, ranks_table.c.height]).execute().fetchall()
        ranks = sorted(ranks, key=lambda x: int(x[1]))  # sort by height
        self.ranks = [r[0] for r in ranks]  # just the ordered ranks
        self.unordered_ranks = set(self.ranks) & set(UNORDERED_RANKS)

        # TODO: can probably remove this check at some point;
        # historically the root node had tax_id == parent_id ==
        # '1', but with a recursive CTE, parent_id must be None for
        # the recursive expression to terminate.
        with self.engine.connect() as con:
            cmd = "select parent_id from {nodes} where rank = 'root'".format(
                nodes=self.nodes)
            result = con.execute(cmd)
            if result.fetchone()[0] is not None:
                raise TaxonIntegrityError('the root node must have parent_id = None')

        self.placeholder = '%s' if self.engine.name == 'postgresql' else '?'

    def _get_table(self, name):
        try:
            val = self.meta.tables[self.prepend_schema(name)]
        except KeyError:
            raise ValueError(
                'Table "{}" is missing from the input database.'.format(name))
        else:
            return val

    def execute(self, statements, exc=IntegrityError, rasie_as=ValueError):
        """Execute ``statements`` in a session, and perform a rollback on
        error. ``exc`` is a single exception object or a tuple of
        objects to be used in the except clause. The error message is
        re-raised as the exception specified by ``raise_as``.

        """

        Session = sessionmaker(bind=self.engine)
        session = Session()

        try:
            for statement in statements:
                session.execute(statement)
        except exc as err:
            session.rollback()
            raise rasie_as(str(err))
        else:
            session.commit()
        finally:
            session.close()

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
        (SELECT new_tax_id FROM {merged}
         WHERE old_tax_id = {x}), {x})
        """.format(x=self.placeholder, merged=self.merged)

        with self.engine.connect() as con:
            result = con.execute(cmd, (tax_id, tax_id))
            return result.fetchone()[0]

    def _get_lineage(self, tax_id, merge_obsolete=True):
        """Return a list of [(rank, tax_id)] describing the lineage of
        tax_id. If ``merge_obsolete`` is True and ``tax_id`` has been
        replaced, use the corresponding value in table merged.

        """

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
          FROM {nodes}
          WHERE tax_id = {}
        UNION ALL
         SELECT p.tax_id, p.parent_id, p.rank
          FROM a JOIN {nodes} p ON a.parent_id = p.tax_id
        )
        SELECT a.rank, a.tax_id FROM a
        JOIN {ranks} using(rank)
        """.format(self.placeholder, nodes=self.nodes, ranks=self.ranks_table)

        # with some versions of sqlite3, an error is raised when no
        # rows are returned; with others, an empty list is returned.
        try:
            with self.engine.connect() as con:
                result = con.execute(cmd, (tax_id,))
                # reorder so that root is first
                lineage = result.fetchall()[::-1]
        except sqlalchemy.exc.ResourceClosedError:
            lineage = []

        if not lineage:
            raise ValueError('tax id "{}" not found'.format(tax_id))

        return lineage

    def prepend_schema(self, name):
        """Prepend schema name to 'name' when a schema is specified

        """
        return '.'.join([self.schema, name]) if self.schema else name

    def _get_lineage_table(self, tax_ids, merge_obsolete=True):
        """Return a list of [(rank, tax_id, tax_name)] describing the lineage
        of tax_id. If ``merge_obsolete`` is True and ``tax_id`` has
        been replaced, use the corresponding value in table merged.

        """

        try:
            with self.engine.connect() as con:
                # insert tax_ids into a temporary table

                temptab = self.prepend_schema(random_name(12))

                cmd = 'CREATE TEMPORARY TABLE "{tab}" (old_tax_id text)'.format(
                    tab=temptab)
                con.execute(cmd)

                log.info('inserting tax_ids into temporary table')
                # TODO: couldn't find an equivalent of "executemany" - does one exist?
                cmd = 'INSERT INTO "{tab}" VALUES ({x})'.format(
                    tab=temptab, x=self.placeholder)
                for tax_id in tax_ids:
                    con.execute(cmd, tax_id)

                log.info('executing recursive CTE')
                cmd = Template("""
                WITH RECURSIVE a AS (
                 SELECT tax_id as tid, 1 AS ord, tax_id, parent_id, rank
                  FROM {{ nodes }}
                  WHERE tax_id in (
                  {% if merge_obsolete %}
                  SELECT COALESCE(m.new_tax_id, "{{ temptab }}".old_tax_id)
                    FROM "{{ temptab }}" LEFT JOIN {{ merged }} m USING(old_tax_id)
                  {% else %}
                  SELECT * from "{{ temptab }}"
                  {% endif %}
                  )
                UNION ALL
                 SELECT a.tid, a.ord + 1, p.tax_id, p.parent_id, p.rank
                  FROM a JOIN {{ nodes }} p ON a.parent_id = p.tax_id
                )
                SELECT a.tid, a.tax_id, a.parent_id, a.rank, tax_name FROM a
                JOIN {{ names }} using(tax_id)
                WHERE names.is_primary
                ORDER BY tid, ord desc
                """).render(
                    temptab=temptab,
                    merge_obsolete=merge_obsolete,
                    merged=self.merged,
                    nodes=self.nodes,
                    names=self.names,
                )

                result = con.execute(cmd)
                rows = result.fetchall()

                con.execute('DROP TABLE "{}"'.format(temptab))
                log.info('returning lineages')
                if not rows:
                    raise ValueError('no tax_ids were found')
                else:
                    returned = {row[0] for row in rows}
                    # TODO: compare set membership, not lengths
                    if len(returned) < len(tax_ids):
                        msg = ('{} tax_ids were provided '
                               'but only {} were returned').format(
                                   len(tax_ids), len(returned))
                        log.error('Input tax_ids not represented in output:')
                        log.error(sorted(set(tax_ids) - returned))
                        raise ValueError(msg)

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
        """Public method for returning a lineage; includes tax_name and rank

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

    def add_source(self, source_name, description=None):
        """Adds a row to table "source" if "name" does not
        exist. Returns (source_id, True) if a new row is created,
        (source_id, False) otherwise.

        """

        # TODO: shoud be able to do this inside a transaction

        if not source_name:
            raise ValueError('"source_name" may not be None or an empty string')

        sel = select([self.source.c.id], self.source.c.name == source_name).execute()
        result = sel.fetchone()
        if result:
            return result[0], False
        else:
            ins = self.source.insert().execute(
                name=source_name, description=description)
            return ins.inserted_primary_key[0], True

    def get_source(self, source_id=None, source_name=None):
        """Returns a dict with keys ['id', 'name', 'description'] or None if
        no match. The ``id`` field is guaranteed to be an int that
        exists in table source. Requires exactly one of ``source_id``
        or ``source_name``. A new source corresponding to
        ``source_name`` is created if necessary.

        """

        if not (bool(source_id) ^ bool(source_name)):
            raise ValueError('exactly one of source_id or source_name is required')

        if source_id:
            try:
                source_id = int(source_id)
            except (ValueError, TypeError):
                raise ValueError(
                    'source_id must be an int or a string representing one')

            sel = select([self.source], self.source.c.id == source_id).execute()
        else:
            sel = select([self.source], self.source.c.name == source_name).execute()

        result = sel.fetchone()
        if not result:
            raise ValueError(
                'there is no source with id {} or name {}'.format(
                    source_id, source_name))

        return dict(list(zip(list(sel.keys()), result)))

    def verify_rank_integrity(self, tax_id, rank, parent_id, children):
        """Confirm that for each node the parent ranks and children ranks are
        coherent

        """
        def _lower(n1, n2):
            return self.ranks.index(n1) < self.ranks.index(n2)

        if rank not in self.ranks:
            raise TaxonIntegrityError('rank "{}" is undefined'.format(rank))

        parent_rank = self.rank(parent_id)
        # undefined ranks can be placed anywhere in a lineage
        if not _lower(rank, parent_rank) and rank not in self.unordered_ranks:
            msg = ('New node "{}", rank "{}" has same or '
                   'higher rank than parent node "{}", rank "{}"')
            msg = msg.format(tax_id, rank, parent_id, parent_rank)
            raise TaxonIntegrityError(msg)

        for child in children:
            if not _lower(self.rank(child), rank):
                msg = 'Child node {} has same or lower rank as new node {}'
                msg = msg.format(tax_id, child)
                raise TaxonIntegrityError(msg)
        return True

    def has_node(self, tax_id):
        result = select([self.nodes], self.nodes.c.tax_id == tax_id).execute()
        return bool(result.fetchone())

    def add_node(self, tax_id, parent_id, rank, names, source_name, children=None,
                 is_valid=True, execute=True, **ignored):
        """Add a node to the taxonomy.

        ``source_name`` is added to table "source" if necessary.

        """

        if ignored:
            log.info('some arguments were ignored: {} '.format(str(ignored)))

        children = children or []
        self.verify_rank_integrity(tax_id, rank, parent_id, children)
        source_id, __ = self.add_source(source_name)
        assert isinstance(is_valid, bool)

        statements = []

        # add node
        statements.append(
            self.nodes.insert().values(
                tax_id=tax_id,
                parent_id=parent_id,
                is_valid=is_valid,
                rank=rank,
                source_id=source_id))

        # add names. Since this is a new node, at least one name must
        # be provided; if only one is provided, it is the primary
        # name. If more than one is primary, an error will be raised
        # from add_names()
        if len(names) == 1:
            names[0]['is_primary'] = True
        else:
            primary_names = [n['tax_name'] for n in names if n.get('is_primary')]
            if len(primary_names) != 1:
                raise ValueError(
                    '`is_primary` must be True for exactly one name in `names`')

        for namedict in names:
            namedict['source_id'] = source_id
            if 'source_name' in namedict:
                del namedict['source_name']

        statements.extend(self.add_names(tax_id, names, execute=False))

        # add children and update source_id
        for child in children:
            statements.append(self.nodes.update(
                whereclause=self.nodes.c.tax_id == child,
                values={'parent_id': tax_id, 'source_id': source_id}))

        if execute:
            self.execute(statements)
        else:
            return statements

    def update_node(self, tax_id, source_name, parent_id=None, rank=None, names=None,
                    children=None, is_valid=None, execute=True, **ignored):

        children = children or []
        source_id, __ = self.add_source(source_name)

        statements = []

        result = select([self.nodes], self.nodes.c.tax_id == tax_id).execute()
        current = dict(list(zip(list(result.keys()), result.fetchone())))

        # update node
        values = dict(source_id=source_id)
        if parent_id is None:
            parent_id = current['parent_id']
        else:
            values['parent_id'] = parent_id

        if rank is not None:
            self.verify_rank_integrity(tax_id, rank, parent_id, children)
            values['rank'] = rank

        if is_valid is not None:
            assert isinstance(is_valid, bool)
            values['is_valid'] = is_valid

        statements.append(self.nodes.update(
            whereclause=self.nodes.c.tax_id == tax_id, values=values))

        # add names if any are provided; names are assumed to be new;
        # there is no checking for primary names.
        if names:
            for namedict in names:
                namedict['source_id'] = source_id
                if 'source_name' in namedict:
                    del namedict['source_name']

            statements.extend(self.add_names(tax_id, names, execute=False))

        # add children
        for child in children:
            statements.append(self.nodes.update(
                whereclause=self.nodes.c.tax_id == child,
                values={'parent_id': tax_id}))

        if execute:
            self.execute(statements)
        else:
            return statements

    def add_name(self, tax_id, tax_name, source_name=None, source_id=None,
                 name_class='synonym', is_primary=False, is_classified=None,
                 execute=True, **ignored):

        """Add a record to the names table corresponding to
        ``tax_id``. Arguments are as follows:

        - tax_id (string, required)
        - tax_name (string, required)

        *one* of the following are required:

        - source_id (int or string coercable to int)
        - source_name (string)

        ``source_id`` or ``source_name`` must identify an existing
        record in table "source".

        The following are optional:

        - name_class (string, default 'synonym')
        - is_primary (bool, see below)
        - is_classified (bool or None, default None)

        ``is_primary`` is optional and defaults to True if only one
        name is provided; otherwise is_primary must be True for
        exactly one name (and is optional in others).

        """

        assert isinstance(is_primary, bool)
        assert is_classified in {None, True, False}
        if ignored:
            log.info('some arguments were ignored: {} '.format(str(ignored)))

        source_id = self.get_source(source_id, source_name)['id']

        statements = []

        if is_primary:
            statements.append(self.names.update(
                whereclause=self.names.c.tax_id == tax_id,
                values={'is_primary': False}))

        statements.append(self.names.insert().values(
            tax_id=tax_id,
            tax_name=tax_name,
            source_id=source_id,
            is_primary=is_primary,
            name_class=name_class,
            is_classified=is_classified))

        if execute:
            self.execute(statements)
        else:
            return statements

    def add_names(self, tax_id, names, execute=True):
        """Associate one or more names with ``tax_id``.

        ``names`` is a list of one or more dicts, with keys
        corresponding to the signature of ``self.add_name()``
        (excluding ``execute``).

        """

        primary_names = [n['tax_name'] for n in names if n.get('is_primary')]
        if len(primary_names) > 1:
            raise ValueError(
                '`is_primary` may be True for no more than one name in `names`')

        statements = []

        for namevals in names:
            if 'tax_id' in namevals:
                del namevals['tax_id']
            statements.extend(
                self.add_name(tax_id=tax_id, execute=False, **namevals))

        if execute:
            self.execute(statements)
        else:
            return statements

    def sibling_of(self, tax_id):
        """Return None or a tax_id of a sibling of *tax_id*.

        If ``tax_id`` is None, then always returns None. Otherwise,
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
        lineage = self.lineage(node)
        return ancestor in list(lineage.values())

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
        # TODO: replace with recursive CTE
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
