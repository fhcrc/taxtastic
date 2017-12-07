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
#    along with taxtastic.  If not, see <http://www.gnu.org/licenses/>
import csv


class Tree(object):
    """Tree for describing taxonomies."""

    def __init__(self, key, **nodedata):
        self.key = key
        self.data = nodedata
        self.parent = None
        self.children = []
        self.descendents = {key: self}

    def __repr__(self, n=0):
        return "  " * n + "Tree(%s" % self.key + "".join(', %s=%s' % (k, v) for k, v in self.data.items()) + ")" + \
            ("" if len(self.children) == 0 else "(\n" +
             ",\n".join(c.__repr__(n + 1) for c in self.children) + ")") + ''

    def __call__(self, *children):
        for c in children:
            c.parent = self
            self.children.append(c)
            # When creating a tree a priori, the children are fully created
            # before the parents are, so we have to add all
            # descendents of the children to the parents when we
            # finally get to them.
            self.descendents.update(c.descendents)
            # On the other hand, when adding children to an existing
            # tree, we have to propogate them up all the descendents
            # dictionaries of the parents.
            q = self
            while True:
                q.descendents[c.key] = c
                if q.isroot():
                    break
                else:
                    q = q.parent
        return self

    def __getattribute__(self, name):
        if name == 'children':
            return object.__getattribute__(self, 'children')
        elif name in object.__getattribute__(self, 'data'):
            return object.__getattribute__(self, 'data')[name]
        else:
            return object.__getattribute__(self, name)

    def isroot(self):
        return self.parent == self or self.parent is None

    def lonelynodes(self):
        return [x for x in self.descendents.values()
                if x.parent is not None and len(x.parent.children) == 1]


def taxtable_to_tree(handle):
    """Read a CSV taxonomy from *handle* into a Tree."""
    c = csv.reader(handle, quoting=csv.QUOTE_NONNUMERIC)
    header = next(c)
    rootdict = dict(list(zip(header, next(c))))
    t = Tree(rootdict['tax_id'], rank=rootdict[
             'rank'], tax_name=rootdict['tax_name'])
    for l in c:
        d = dict(list(zip(header, l)))
        target = t.descendents[d['parent_id']]
        target(Tree(d['tax_id'], rank=d['rank'], tax_name=d['tax_name']))
    return t


def lonely_company(taxonomy, tax_ids):
    """Return a set of species tax_ids which will makes those in *tax_ids* not lonely.

    The returned species will probably themselves be lonely.
    """
    return [taxonomy.species_below(taxonomy.sibling_of(t)) for t in tax_ids]


def solid_company(taxonomy, tax_ids):
    """Return a set of non-lonely species tax_ids that will make those in *tax_ids* not lonely."""
    res = []
    for t in tax_ids:
        res.extend(taxonomy.nary_subtree(taxonomy.sibling_of(t), 2) or [])
    return res
