import algotax

_leaf_numbers = {}

class Leaf(object):
    is_leaf = True
    def __init__(self, color):
        self.color = color
        self.parent = None
        self.number = _leaf_numbers[self] = len(_leaf_numbers)

    def __repr__(self):
        return '<%r#%d>' % (self.color, self.number)

class Edge(object):
    def __init__(self, top, bottom):
        self.top = top
        self.bottom = bottom
        self.cut_colors = set()

    def __repr__(self):
        colors = ''
        if self.cut_colors:
            colors = '[%s]' % '; '.join(str(c) for c in self.cut_colors)
        return '-%s> %r' % (colors, self.bottom)

class Node(object):
    is_leaf = False
    def __init__(self, lft, rgt):
        self.lft = lft.parent = Edge(self, lft)
        self.rgt = rgt.parent = Edge(self, rgt)
        self.parent = None

    @property
    def children(self):
        return self.lft, self.rgt

    def __repr__(self):
        return '(%r, %r)' % (self.lft, self.rgt)

def build_tree(tree):
    def _aux(cur):
        if len(cur) == 1:
            return cur[0]
        lft, rgt = map(_aux, cur)
        return Node(lft, rgt)
    return _aux(tree)

def intersection(it):
    ret = None
    for x in it:
        if ret is None:
            ret = set(x)
        else:
            ret.intersection_update(x)
    return ret or set()

def color_edges(root):
    stack = [('down', root, None)]
    while stack:
        phase, cur, color = stack.pop()
        if phase == 'down':
            if cur.is_leaf:
                stack.append(('up', cur.parent, cur.color))
            else:
                stack.extend(('down', e.bottom, None) for e in cur.children)
        elif phase == 'up':
            if cur is None or color in cur.cut_colors:
                continue
            cur.cut_colors.add(color)
            stack.append(('up', cur.top.parent, color))

    # This isn't actually generalized to >2 children.
    # It's just shorter this way.
    stack = [(root, set())]
    while stack:
        node, okayed = stack.pop()
        if node.is_leaf:
            continue
        okayed = intersection(e.cut_colors for e in node.children) | okayed
        for e in node.children:
            e_ = e.cut_colors & okayed
            if e_ != e.cut_colors:
                stack.append((e.bottom, okayed))
                e.cut_colors = e_

def main():
    root = build_tree([
        [[Leaf(1)], [Leaf(2)]], [[Leaf(1)], [[Leaf(2)], [Leaf(2)]]]
    ])
    print root
    color_edges(root)
    print root
    print algotax.walk(root)

if __name__ == '__main__':
    main()
