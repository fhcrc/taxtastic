from Bio import Phylo
from StringIO import StringIO
import unittest
from taxtastic import algotax

class ColoredTreeTestMixin(object):
    @classmethod
    def setup_class(cls):
        cls.parsed_tree = Phylo.read(StringIO(cls.tree), 'newick')
        cls.colors = {n: n.name for n in cls.parsed_tree.get_terminals()}
        cls.metadata = algotax.color_clades(cls.parsed_tree, cls.colors)

    @classmethod
    def setUpClass(cls):
        super(ColoredTreeTestMixin, cls).setUpClass()
        cls.setup_class()

class CladeColorTestMixin(ColoredTreeTestMixin):
    def test_parents(self):
        self.assertEqual(
            len(self.metadata.parents),
            sum((1 for _ in self.parsed_tree.find_clades())))
        for child, parent in self.metadata.parents.iteritems():
            if parent is None:
                self.assertIs(child, self.parsed_tree.root)
            else:
                self.assertIn(child, parent.clades)

    def test_cut_colors(self):
        numbering = dict(enumerate(
            self.parsed_tree.find_clades(order='postorder')))
        rev_numbering = {v: k for k, v in numbering.iteritems()}
        for color, nodes in self.cut_colors.iteritems():
            if color is None:
                for node in nodes:
                    cut_colors = self.metadata.cut_colors[numbering[node]]
                    self.assert_(not cut_colors)
                continue

            for node in numbering.itervalues():
                if node is self.parsed_tree.root:
                    continue
                self.assertEqual(
                    rev_numbering[node] in nodes,
                    color in self.metadata.cut_colors[node])

class CladeColorTest1(CladeColorTestMixin, unittest.TestCase):
    tree = '((A,A),(B,B))'
    cut_colors = {
        'A': {0, 1},
        'B': {3, 4},
        None: {2, 5},
    }

class CladeColorTest2(CladeColorTestMixin, unittest.TestCase):
    tree = '((A,B),((A,B),A))'
    cut_colors = {
        'A': {0, 2, 7, 5, 3, 6},
        'B': {1, 2, 7, 5, 4},
    }

class CladeColorTest3(CladeColorTestMixin, unittest.TestCase):
    tree = '(((A,B),(C,D)),(A,B),(C,A))'
    cut_colors = {
        'A': {0, 2, 6, 9, 7, 12, 11},
        'B': {1, 2, 6, 9, 8},
        'C': {3, 5, 6, 12, 10},
        'D': set(),
        None: {4},
    }

class AlgotaxWalkTestMixin(ColoredTreeTestMixin):
    @classmethod
    def setup_class(cls):
        super(AlgotaxWalkTestMixin, cls).setup_class()
        cls.nodeset = algotax.walk(cls.parsed_tree.root, cls.metadata)

    setUpClass = setup_class

    def test_walk(self):
        self.assertEqual(len(self.nodeset), self.convex_tree_size)

class AlgotaxWalkTest1(AlgotaxWalkTestMixin, unittest.TestCase):
    tree = '((A,A),(B,B))'
    convex_tree_size = 4

class AlgotaxWalkTest2(AlgotaxWalkTestMixin, unittest.TestCase):
    tree = '((A,B),(A,B))'
    convex_tree_size = 3

class AlgotaxWalkTest3(AlgotaxWalkTestMixin, unittest.TestCase):
    tree = '(((A,B),(A,B)),(C,C))'
    convex_tree_size = 5

class AlgotaxWalkTest4(AlgotaxWalkTestMixin, unittest.TestCase):
    tree = '(((A,B),B),(A,A))'
    convex_tree_size = 4

class RerootingTestMixin(object):
    @classmethod
    def setup_class(cls):
        cls.parsed_tree = Phylo.read(StringIO(cls.tree), 'newick')
        cls.mrcas = {n: int(n.branch_length)
            for n in cls.parsed_tree.get_terminals()}
        cls.new_root = algotax.reroot(cls.parsed_tree.root, cls.mrcas.get)

    @classmethod
    def setUpClass(cls):
        super(RerootingTestMixin, cls).setUpClass()
        cls.setup_class()

    def test_rerooting(self):
        root_number = next(e
            for e, n in enumerate(
                self.parsed_tree.find_clades(order='postorder'))
            if n is self.new_root)
        self.assertEqual(root_number, self.root_number)

class RerootingTest1(RerootingTestMixin, unittest.TestCase):
    tree = '(0,0)'
    root_number = 2

class RerootingTest2(RerootingTestMixin, unittest.TestCase):
    tree = '(0,(2,2)0)'
    root_number = 3

class RerootingTest3(RerootingTestMixin, unittest.TestCase):
    tree = '(6,(2,((7,7)3,(7,7)3)0),6)'
    root_number = 8

class RerootingTest4(RerootingTestMixin, unittest.TestCase):
    tree = '((((6,7)4,5)2,3)0,1)'
    root_number = 0
