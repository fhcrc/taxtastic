import csv

class Tree(object):
    def __init__(self, key, **nodedata):
        self.key = key
        self.data = nodedata
        self.parent = None
        self.children = []
        self.descendents = {key: self}
    def __repr__(self, n=0):
        return "  "*n + "Tree(%s" % self.key + "".join(', %s=%s' % (k,v) for k,v in self.data.iteritems()) + ")" + \
            ("" if len(self.children) == 0 else "(\n" + ",\n".join(c.__repr__(n+1) for c in self.children) + ")") + ''
    def __call__(self, *children):
        for c in children:
            c.parent = self
            self.children.append(c)
            self.descendents.update(c.descendents)
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
        return [x for x in self.descendents.itervalues()
                if x.parent is not None and len(x.parent.children) == 1]


def taxtable_to_tree(handle):
    c = csv.reader(handle, quoting=csv.QUOTE_NONNUMERIC)
    header = c.next()
    rootdict = dict(zip(header, c.next()))
    t = Tree(rootdict['tax_id'], rank=rootdict['rank'], tax_name=rootdict['tax_name'])
    for l in c:
        d = dict(zip(header, l))
        target = t.descendents[d['parent_id']]
        target(Tree(d['tax_id'], rank=d['rank'], tax_name=d['tax_name']))
    return t

    
