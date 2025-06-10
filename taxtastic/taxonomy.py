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

import sqlalchemy as sa
from sqlalchemy import MetaData, and_, or_
from sqlalchemy.sql import select
from sqlalchemy.exc import IntegrityError
from sqlalchemy.orm import Session

from taxtastic.ncbi import UnknownRankError, UNORDERED_RANKS
from taxtastic.utils import random_name

log = logging.getLogger(__name__)


class TaxonIntegrityError(Exception):
    '''
    Raised when something in the Taxonomy is not structured correctly
    '''


class Taxonomy(object):

    def __init__(self, engine, schema=None):
        """The Taxonomy class defines an object providing an interface
        to the taxonomy database.

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
        self.meta.reflect(bind=self.engine)
        self.schema = schema

        self.nodes = self._get_table('nodes')
        self.names = self._get_table('names')
        self.source = self._get_table('source')
        self.merged = self._get_table('merged')
        ranks_table = self._get_table('ranks')
        self.ranks_table = ranks_table

        # intended for constructing text queries using the pattern
        # cmd = 'select * from {<table>}'.format(**self.tablenames)
        # to ensure that the schema is prepended to the table name when defined
        self.tables = {name: self.prepend_schema(name) for name in [
            'nodes', 'names', 'source', 'merged', 'ranks']}

        ranks = self.fetchall(
            select(self.ranks_table.c.rank).order_by(ranks_table.c.height))
        self.ranks = [r[0] for r in ranks]

        self.unordered_ranks = set(self.ranks) & set(UNORDERED_RANKS)

        # parent_id must be None for the recursive CTE for calculating
        # lineages to terminate.
        parent_node = self.fetchone(
            select(self.nodes).filter_by(rank='root'))

        if parent_node.parent_id is not None:
            raise TaxonIntegrityError(
                'the root node must have parent_id = None')

        self.placeholder = '%s' if self.engine.name == 'postgresql' else '?'

    def _get_table(self, name):
        try:
            val = self.meta.tables[self.prepend_schema(name)]
        except KeyError as ex:
            raise ValueError(
                f'Table "{name}" is missing from the input database.') from ex
        else:
            return val

    def execute(self, statements, exc=IntegrityError, raise_as=ValueError,
                errormsg=''):
        """Execute list of ``statements`` in a transaction, and
        perform a rollback on error. ``exc`` is a single exception
        object or a tuple of objects to be used in the except clause.
        The error message is re-raised as the exception specified by
        ``raise_as``.

        """

        try:
            with self.engine.begin() as conn:
                for stmt in statements:
                    conn.execute(stmt)
        except exc as ex:
            raise raise_as(errormsg) from ex

    def _node(self, tax_id):
        """
        Returns parent_id, rank

        FIXME: expand return rank to include custom 'below' ranks built when
               get_lineage is caled
        """

        output = self.fetchone(
            select(self.nodes.c.parent_id, self.nodes.c.rank)
            .filter_by(tax_id=tax_id))

        if not output:
            raise ValueError(f'value "{tax_id}" not found in nodes.tax_id')
        else:
            return output

    def id_from_names(self, tax_names):
        """Return tax_ids corresponding to tax_names.

        """
        names = self.names
        s = select(names.c.tax_name, names.c.tax_id).\
            where(names.c.tax_name.in_(tax_names)).\
            distinct()
        return self.fetchall(s)

    def primary_from_id(self, tax_id):
        """Returns primary taxonomic name associated with tax_id

        """

        output = self.fetchone(
            select(self.names.c.tax_name)
            .where(and_(self.names.c.tax_id == tax_id,
                        self.names.c.is_primary)))

        if output:
            return output[0]
        else:
            raise ValueError(f'"{tax_id}" not found in names.tax_id')

    def primary_from_ids(self, tax_ids):
        names = self.names
        s = select(names.c.tax_id, names.c.tax_name).\
            where(and_(
                names.c.tax_id.in_(tax_ids),
                names.c.is_primary)).\
            distinct()
        return self.fetchall(s)

    def primary_from_name(self, tax_name):
        """
        Return tax_id and primary tax_name corresponding to tax_name.
        """

        res = self.fetchone(
            select(self.names.c.tax_id, self.names.c.is_primary)
            .filter_by(tax_name=tax_name))

        if res:
            tax_id, is_primary = res
        else:
            msg = '"{}" not found in names.tax_names'.format(tax_name)
            raise ValueError(msg)

        if not is_primary:
            tax_name = self.primary_from_id(self.names.c.tax_id)

        return tax_id, tax_name, bool(is_primary)

    def _get_merged(self, tax_id):
        """Returns tax_id into which `tax_id` has been merged or
        `tax_id` if not obsolete.

        """

        cmd = sa.text("""
        SELECT COALESCE(
        (SELECT new_tax_id FROM {merged}
         WHERE old_tax_id = :tax_id), :tax_id)
        """.format(merged=self._get_table('merged')))

        return self.fetchone(cmd, tax_id=tax_id)[0]

    def _get_lineage(self, tax_id, merge_obsolete=True):
        """Return a list of [(rank, tax_id)] describing the lineage of
        tax_id from root to tip. If ``merge_obsolete`` is True and
        ``tax_id`` has been replaced, use the corresponding value in
        table merged.

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

        nodes = self._get_table('nodes')
        ranks = self._get_table('ranks')

        cmd = sa.text(f"""
        WITH RECURSIVE a AS (
         SELECT tax_id, parent_id, rank
          FROM {nodes}
          WHERE tax_id = :tax_id
        UNION ALL
         SELECT p.tax_id, p.parent_id, p.rank
          FROM a JOIN {nodes} p ON a.parent_id = p.tax_id
        )
        SELECT a.rank, a.tax_id FROM a
        JOIN {ranks} using(rank)
        """)

        # with some versions of sqlite3, an error is raised when no
        # rows are returned; with others, an empty list is returned.
        try:
            lineage = self.fetchall(cmd, tax_id=tax_id)[::-1]
        except sa.exc.ResourceClosedError:
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

                cmd = f'CREATE TEMPORARY TABLE "{temptab}" (old_tax_id text)'
                con.execute(sa.text(cmd))

                log.info('inserting tax_ids into temporary table')
                cmd = sa.text(f'INSERT INTO "{temptab}" VALUES (:tax_id)')
                con.execute(cmd, [{'tax_id': tax_id} for tax_id in tax_ids])

                log.info('executing recursive CTE')
                cmd = sa.text(Template("""
                WITH RECURSIVE a AS (
                 SELECT tax_id as tid, 1 AS ord, tax_id, parent_id, rank
                  FROM "{{ nodes }}"
                  WHERE tax_id in (
                  {% if merge_obsolete %}
                  SELECT COALESCE(m.new_tax_id, "{{ temptab }}".old_tax_id)
                    FROM "{{ temptab }}"
                    LEFT JOIN {{ merged }} m USING(old_tax_id)
                  {% else %}
                  SELECT * from "{{ temptab }}"
                  {% endif %})
                UNION ALL
                 SELECT a.tid, a.ord + 1, p.tax_id, p.parent_id, p.rank
                  FROM a JOIN "{{ nodes }}" p ON a.parent_id = p.tax_id
                )
                SELECT a.tid, a.tax_id, a.parent_id, a.rank, tax_name FROM a
                JOIN "{{ names }}" using(tax_id)
                WHERE names.is_primary
                ORDER BY tid, ord desc
                """).render(
                    temptab=temptab,
                    merge_obsolete=merge_obsolete,
                    **self.tables))

                result = con.execute(cmd)
                rows = result.fetchall()

                con.execute(sa.text(f'DROP TABLE "{temptab}"'))

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

        except sa.exc.ResourceClosedError:
            raise ValueError('tax id "{}" not found'.format(tax_id))

    def is_below(self, lower, upper):
        return lower in self.ranks_below(upper)

    def ranks_below(self, rank, depth=None):
        below = []
        try:
            index = self.ranks.index(rank)
            depth = 0 if depth is None else index - depth
            below = self.ranks[depth:index]
        except ValueError as err:
            log.error(err)
        return below

    def synonyms(self, tax_id=None, tax_name=None):
        if not bool(tax_id) ^ bool(tax_name):
            raise ValueError(
                'Exactly one of tax_id and tax_name may be provided.')

        names = self.names

        if tax_name:
            res = self.fetchone(
                select(names.c.tax_id)
                .filter_by(tax_name=tax_name))

            if res:
                tax_id = res[0]
            else:
                msg = f'"{tax_name}" not found in names.tax_names'
                raise ValueError(msg)

        output = self.fetchall(
            select(names.c.tax_name, names.c.is_primary)
            .filter_by(tax_id=tax_id))

        if not output:
            raise ValueError(f'"{tax_name}" not found in names.tax_id')

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

    def fetchone(self, statement, **params):
        """Return first result of select statement 'statement'. If
        statement is a text statement, variable substitutions can be
        provided as keyword arguments.

        """
        with Session(self.engine) as session:
            return session.execute(statement, params).fetchone()

    def fetchall(self, statement, **params):
        """Return all results of select statement 'statement'. If
        statement is a text statement, variable substitutions can be
        provided as keyword arguments.

        """
        with Session(self.engine) as session:
            return session.execute(statement, params).fetchall()

    def add_source(self, source_name, description=None):
        """Adds a row to table "source" if "name" does not
        exist. Returns (source_id, True) if a new row is created,
        (source_id, False) otherwise.

        """

        if not source_name:
            raise ValueError('"source_name" may not be None or an empty string')

        result = self.fetchone(
            select(self.source.c.id).filter_by(name=source_name))

        if result:
            return (result[0], False)
        else:
            with self.engine.begin() as conn:
                stmt = sa.insert(self.source).values(
                    name=source_name, description=description)
                result = conn.execute(stmt)
                return (result.inserted_primary_key[0], True)

    def get_source(self, source_id=None, source_name=None):
        """Returns a dict with keys ('id', 'name', 'description'). The
        ``id`` field is guaranteed to be an int that exists in table
        source. Requires exactly one of ``source_id`` or
        ``source_name``.

        """

        if not (bool(source_id) ^ bool(source_name)):
            raise ValueError(
                'exactly one of source_id and source_name is required')

        if source_id:
            try:
                condition = {'id': int(source_id)}
            except (ValueError, TypeError) as ex:
                raise ValueError(f'"{source_id}" is not an integer') from ex
        elif source_name:
            condition = {'name': source_name}

        result = self.fetchone(select(self.source).filter_by(**condition))

        if not result:
            raise ValueError(
                'there is no source with id {} or name {}'.format(
                    source_id, source_name))

        return result._asdict()

    def verify_rank_integrity(self, tax_id, rank, parent_id, children):
        """Confirm that for each node the parent ranks and children ranks are
        coherent

        """
        def _lower(n1, n2):
            return self.ranks.index(n1) < self.ranks.index(n2)

        if rank not in self.ranks:
            raise UnknownRankError('"{}" is undefined'.format(rank))

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
        result = self.fetchone(
            select(self.nodes)
            .filter_by(tax_id=tax_id))
        return bool(result)

    def unknowns(self, tax_ids):
        result = self.fetchall(
            select(self.nodes.c.tax_id)
            .filter(self.nodes.c.tax_id.in_(tax_ids)))
        result = set(r for r, in result)
        return [i for i in tax_ids if i not in result]

    def add_node(self, tax_id, parent_id, rank, names, source_name,
                 children=None, is_valid=True, execute=True, **ignored):
        """Add a node to the taxonomy.

        ``source_name`` is added to table "source" if necessary.
        """

        if ignored:
            log.info(f'some arguments were ignored: {ignored} ')

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
            statements.append(
                sa.update(self.nodes)
                .where(self.nodes.c.tax_id == child)
                .values(parent_id=tax_id, source_id=source_id))

        if execute:
            self.execute(statements)
        else:
            return statements

    def update_node(self, tax_id, source_name, parent_id=None, rank=None,
                    names=None, children=None, is_valid=None, execute=True,
                    **ignored):

        children = children or []
        source_id, __ = self.add_source(source_name)

        statements = []

        result = self.fetchone(
            select(self.nodes)
            .where(self.nodes.c.tax_id == tax_id))
        current = result._asdict()

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

        statements.append(
            sa.update(self.nodes)
            .where(self.nodes.c.tax_id == tax_id)
            .values(values))

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
            statements.append(
                sa.update(self.nodes)
                .where(self.nodes.c.tax_id == child)
                .values(parent_id=tax_id))

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
            log.info(f'some arguments were ignored: {ignored} ')

        source_id = self.get_source(source_id, source_name)['id']

        statements = []

        if is_primary:
            statements.append(
                sa.update(self.names)
                .where(self.names.c.tax_id == tax_id)
                .values(is_primary=False))

        statements.append(
            sa.insert(self.names)
            .values(tax_id=tax_id,
                    tax_name=tax_name,
                    source_id=source_id,
                    is_primary=is_primary,
                    name_class=name_class,
                    is_classified=is_classified))

        if execute:
            self.execute(
                statements, errormsg=f'{tax_id}, {tax_name} already exists')
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
        parent_id, rank = self._node(tax_id)

        output = self.fetchone(
            select(self.nodes.c.tax_id)
            .where(and_(self.nodes.c.parent_id == parent_id,
                        self.nodes.c.tax_id != tax_id,
                        self.nodes.c.rank == rank)))

        if output:
            return output[0]
        else:
            log.info(
                f'No sibling of tax_id {tax_id} '
                f'with rank {rank} found in taxonomy')
            return None

    def is_ancestor_of(self, node, ancestor):
        if node is None or ancestor is None:
            return False
        lineage = self.lineage(node)
        return ancestor in list(lineage.values())

    def rank(self, tax_id):
        return self._node(tax_id)[1]

    def tax_ids(self):
        '''
        Return all tax_ids in node table
        '''
        fetch = self.fetchall(select(self.nodes.c.tax_id))
        ids = [t[0] for t in fetch]
        return ids

    def child_of(self, tax_id):
        """Return None or a tax id of a child of *tax_id*.

        If *tax_id* is None, then always returns None. Otherwise
        returns a child if one exists, else None. The child must have
        a proper rank below that of tax_id (i.e., genus, species, but
        not no_rank or below_below_kingdom).
        """
        __, rank = self._node(tax_id)

        output = self.fetchone(
            select(self.nodes.c.tax_id)
            .where(and_(self.nodes.c.parent_id == tax_id,
                        or_(*[self.nodes.c.rank == r
                              for r in self.ranks_below(rank)]))))

        if output:
            r = output[0]
            assert self.is_ancestor_of(r, tax_id)
            return r
        else:
            log.warning(f'No children of tax_id {tax_id} with '
                        f'rank below {rank} found in database')
            return None

    def children_of(self, tax_id, n):

        # TODO: replace with recursive CTE?
        __, rank = self._node(tax_id)
        output = self.fetchall(
            select(self.nodes.c.tax_id)
            .where(and_(self.nodes.c.parent_id == tax_id,
                        or_(*[self.nodes.c.rank == r
                              for r in self.ranks_below(rank)])))
            .limit(n))

        if output:
            r = [x[0] for x in output]
            for x in r:
                assert self.is_ancestor_of(x, tax_id)
            return r
        else:
            return []

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
        _, rank = self._node(tax_id)
        if rank == 'species':
            return [tax_id]
        else:
            children = self.children_of(tax_id, 2)
            species_taxids = []
            for t in children:
                species_taxids.extend(self.nary_subtree(t, n))
            return species_taxids

    def species_below(self, tax_id):
        try:
            _, rank = self._node(tax_id)
        except ValueError:
            return None
        if rank == 'species':
            return tax_id
        else:
            c = self.child_of(tax_id)
            newc = self.species_below(c)
            assert self.is_ancestor_of(newc, tax_id)
            return newc

    def descendants_of(self, tax_ids):
        """Return list of all tax_ids under *tax_id*"""
        tax_ids = ','.join("'{}'".format(t) for t in tax_ids)
        cmd = sa.text("""
        WITH RECURSIVE descendants AS (
         SELECT tax_id
         FROM nodes
         WHERE tax_id in ({})
         UNION ALL
         SELECT
         n.tax_id
         FROM nodes n
         JOIN descendants d ON d.tax_id = n.parent_id
        ) SELECT DISTINCT tax_id
        FROM descendants
        JOIN names using(tax_id)
        WHERE is_primary;
        """.format(tax_ids))
        with self.engine.connect() as con:
            return [row[0] for row in con.execute(cmd).fetchall()]

    def is_valid(self, tax_ids=None, no_rank=True):
        """Return all classified tax_ids"""
        nodes = self.nodes
        s = select(nodes.c.tax_id).where(nodes.c.is_valid)
        if tax_ids:
            s = s.where(nodes.c.tax_id.in_(set(tax_ids)))
        if not no_rank:
            s = s.where(nodes.c.rank == 'no_rank')
        return [r[0] for r in self.fetchall(s)]
