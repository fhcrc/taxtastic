import os.path
import unittest

try:
    # python2
    from cStringIO import StringIO
except ImportError:
    # python3
    from io import StringIO

from taxtastic.taxtable import TaxNode
from .config import data_path

DN = os.path.dirname(__file__)


class TaxNodeTestCase(unittest.TestCase):

    def setUp(self):
        with open(data_path('simple_taxtable.csv')) as fp:
            self.root = TaxNode.from_taxtable(fp)

    def test_index(self):
        self.assertEqual(self.root, self.root.get_node('1'))
        self.assertEqual('1', self.root.tax_id)

    def test_iter(self):
        self.assertEqual(356, sum(1 for i in self.root))

    def test_lineage(self):
        node = self.root.get_node('1303')
        lineage = node.lineage()
        self.assertEqual(['1', '131567', '2', '1239', '91061', '186826', '1300', '1301', '1303'],
                         [i.tax_id for i in lineage])

    def test_write_taxtable(self):
        expected = '''"tax_id","parent_id","rank","tax_name","root","cellular_root","domain","phylum","class","order","family","genus","species"
"1","1","root","root","1","","","","","","","",""
"131567","1","cellular_root","cellular organisms","1","131567","","","","","","",""
"2","131567","domain","Bacteria","1","131567","2","","","","","",""
"1239","2","phylum","Firmicutes","1","131567","2","1239","","","","",""
"91061","1239","class","Bacilli","1","131567","2","1239","91061","","","",""
"186826","91061","order","Lactobacillales","1","131567","2","1239","91061","186826","","",""
"1300","186826","family","Streptococcaceae","1","131567","2","1239","91061","186826","1300","",""
"1301","1300","genus","Streptococcus","1","131567","2","1239","91061","186826","1300","1301",""
"1303","1301","species","Streptococcus oralis","1","131567","2","1239","91061","186826","1300","1301","1303"
'''
        node = self.root.get_node('1303')
        s = StringIO()

        node.write_taxtable(s)
        v = s.getvalue()
        self.assertEqual(expected, v)

    def test_prune_unrepresented(self):
        self.root.get_node('1303').sequence_ids.add('sequence1')
        self.root.prune_unrepresented()
        self.assertEqual(set(['1', '131567', '2', '1239', '91061', '186826', '1300', '1301', '1303']),
                         set(self.root.index))

    def test_collapse(self):
        node = self.root.get_node('1300')
        self.root.get_node('1303').sequence_ids.add('seq1')
        self.root.get_node('1301').sequence_ids.add('seq2')
        node.sequence_ids.add('seq3')
        node.collapse(False)
        self.assertEqual(node.sequence_ids, set(['seq1', 'seq2', 'seq3']))
        self.assertEqual(self.root.get_node('1303').sequence_ids, set())
        self.assertEqual(self.root.get_node('1301').sequence_ids, set())

    def test_drop(self):
        tax_id = "1301"
        sequence_ids = ['dsequence1', 'dsequence2']
        to_drop = self.root.get_node(tax_id)
        for i in sequence_ids:
            to_drop.sequence_ids.add(i)
        children = to_drop.children
        parent = to_drop.parent
        to_drop.drop()
        self.assertNotIn(tax_id, self.root.index)
        self.assertNotIn(to_drop, parent.children)
        for child in children:
            self.assertIn(child, parent.children)
            self.assertIn(child.tax_id, self.root.index)
            self.assertEqual(parent, child.parent)
            self.assertEqual(parent.index, child.index)

        self.assertIsNone(to_drop.index)
        self.assertIsNone(to_drop.parent)
        for i in sequence_ids:
            self.assertIn(i, parent.sequence_ids)

    def test_drop_root(self):
        self.assertRaises(ValueError, self.root.drop)
