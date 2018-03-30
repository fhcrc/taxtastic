import itertools

try:
    from pysqlite2 import dbapi2 as sqlite3
except ImportError:
    import sqlite3


class OnUpdate(object):

    def __init__(self, proxied):
        self.proxied = proxied
        self.setter = None

    def __get__(self, inst, cls):
        if inst is None:
            return self
        return getattr(inst, self.proxied)

    def __set__(self, inst, value):
        if self.setter:
            self.setter(inst, value)
        setattr(inst, self.proxied, value)

    def __call__(self, f):
        self.setter = f
        return self


class _IntermediateTaxon(object):

    def __init__(self, tax_id, parent, rank, tax_name):
        self.children = set()
        self.tax_id = tax_id
        self.parent = parent
        self.rank = rank
        self.tax_name = tax_name

    _parent = _adjacent_to = None

    @OnUpdate('_parent')
    def parent(self, p):
        if self.parent is not None:
            self.parent.children.discard(self)
        if p is not None:
            p.children.add(self)

    @OnUpdate('_adjacent_to')
    def adjacent_to(self, n):
        if n is not None:
            self.parent = n.parent

    def iterate_children(self, on_pop=None, including_self=True):
        search_stack = [(None, set([self]))]
        while search_stack:
            if not search_stack[-1][1]:
                parent, _ = search_stack.pop()
                if on_pop:
                    on_pop(parent)
                continue
            node = search_stack[-1][1].pop()
            if node is not self or including_self:
                yield node
            search_stack.append((node, node.children.copy()))


class Taxdb(object):

    def __init__(self, sqlite_db=None):
        if sqlite_db is None:
            sqlite_db = sqlite3.connect(':memory:')
        self.db = sqlite_db

    def __getattr__(self, attr):
        return getattr(self.db, attr)

    def create_tables(self):
        curs = self.db.cursor()

        curs.execute("""
            CREATE TABLE ranks (
              rank TEXT PRIMARY KEY NOT NULL,
              rank_order INT
            )
        """)

        curs.execute("""
            CREATE TABLE taxa (
              tax_id TEXT PRIMARY KEY NOT NULL,
              tax_name TEXT NOT NULL,
              rank TEXT REFERENCES ranks (rank) NOT NULL
            )
        """)

        curs.execute("""
            CREATE TABLE sequences (
              seqname TEXT PRIMARY KEY NOT NULL,
              tax_id TEXT REFERENCES taxa (tax_id) NOT NULL
            )
        """)

        curs.execute("""
            CREATE TABLE hierarchy (
              tax_id TEXT REFERENCES taxa (tax_id) PRIMARY KEY NOT NULL,
              lft INT NOT NULL UNIQUE,
              rgt INT NOT NULL UNIQUE
            )
        """)

        curs.execute("""
            CREATE VIEW parents AS
            SELECT h1.tax_id AS child,
                   h2.tax_id AS parent
            FROM   hierarchy h1
                   JOIN hierarchy h2
                     ON h1.lft BETWEEN h2.lft AND h2.rgt
        """)

    def insert_from_taxtable(self, fieldnames_cb, table):
        curs = self.db.cursor()

        taxon_map = {}
        for row in table:
            parent = taxon_map.get(row['parent_id'])
            taxon = _IntermediateTaxon(
                row['tax_id'], parent, row['rank'], row['tax_name'])
            taxon_map[taxon.tax_id] = taxon

        root = next(iter(taxon_map.values()))
        while root.parent is not None:
            root = root.parent
        counter = itertools.count(1).__next__

        def on_pop(parent):
            if parent is not None:
                parent.rgt = counter()
        for node in root.iterate_children(on_pop=on_pop):
            node.lft = counter()

        fieldnames = fieldnames_cb()
        curs.executemany("INSERT INTO ranks (rank_order, rank) VALUES (?, ?)",
                         enumerate(fieldnames[4:]))
        curs.executemany("INSERT INTO taxa VALUES (?, ?, ?)",
                         ((t.tax_id, t.tax_name, t.rank) for t in taxon_map.values()))
        curs.executemany("INSERT INTO hierarchy VALUES (?, ?, ?)",
                         ((t.tax_id, t.lft, t.rgt) for t in taxon_map.values()))
        self.db.commit()
