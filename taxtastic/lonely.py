class Tree(object):
    def __init__(self, key, **nodedata):
        self.key = key
        self.data = nodedata
        self.parent = None
        self.children = []
        self.descendents = {key: self}
    def __call__(self, *children):
        for c in children:
            c.parent = self
            self.children.append(c)
            q = self
            while q is not None:
                q.descendents[c.key] = c
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

def taxtable_to_tree(handle):
    c = csv.reader(handle, quoting=csv.QUOTE_NONNUMERIC)
    header = c.next()
    rootdict = dict(zip(header, c.next()))
    t = Tree(rootdict['tax_id'], rank=rootdict['rank'], tax_name=rootdict['tax_name'])
    for l in c:
        d = dict(zip(header, l))
        t.descendents[d['parent_id']](Tree(d['tax_id'], rank=d['tank'], tax_name=d['tax_name']))
    return t
    
