from cStringIO import StringIO
import os.path
import unittest

from taxtastic.taxtable import TaxNode
from .config import data_path

DN = os.path.dirname(__file__)


class NodeTestCase(unittest.TestCase):
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
        expected = '''"tax_id","parent_id","rank","tax_name","root","below_root","superkingdom","phylum","class","order","family","genus","species"
"1","1","root","root","1","","","","","","","",""
"131567","1","below_root","cellular organisms","1","131567","","","","","","",""
"2","131567","superkingdom","Bacteria","1","131567","2","","","","","",""
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
        self.assertEquals(expected, v)

