
import logging
import os
import json
import unittest

#import taxtastic
import taxtastic.utils
from . import config


log = logging


outputdir = os.path.abspath(config.outputdir)
datadir = os.path.abspath(config.datadir)

if hasattr(taxtastic.utils, 'read_spreadsheet'):
    class TestReadSpreadsheet(unittest.TestCase):

        def setUp(self):
            self.funcname = '_'.join(self.id().split('.')[-2:])

        def tearDown(self):
            pass

        def test01(self):
            headers, rows = taxtastic.utils.read_spreadsheet(
                os.path.join(datadir,'new_taxa.xls'))
            check = lambda val: isinstance(val, float)
            self.assertTrue(all([check(row['parent_id']) for row in rows]))

        def test02(self):
            headers, rows = taxtastic.utils.read_spreadsheet(
                os.path.join(datadir,'new_taxa.xls'),
                fmts={'tax_id':'%i','parent_id':'%i'}
                )
            check = lambda val: isinstance(val, str) and '.' not in val
            self.assertTrue(all([check(row['parent_id']) for row in rows]))

## TODO: need to test creation of new nodes with csv (as opposed to xls) output
class TestGetNewNodes(unittest.TestCase):

    xlrd_is_installed = hasattr(taxtastic.utils, 'read_spreadsheet')

    def setUp(self):
        self.funcname = '_'.join(self.id().split('.')[-2:])

    def tearDown(self):
        pass

    def check_parent_id(self, rows):
        check = lambda val: isinstance(val, str) and '.' not in val
        self.assertTrue(all([check(d['parent_id']) for d in rows]))

    def check_children(self, rows):
        for d in rows:
            if 'children' in d:
                val = d['children']
                self.assertTrue(val and isinstance(val, list))

    def test01(self):
        if self.xlrd_is_installed:
            rows = taxtastic.utils.get_new_nodes(os.path.join(datadir,'new_taxa.xls'))
            check = lambda val: isinstance(val, str) and '.' not in val
            self.assertTrue(all([check(row['parent_id']) for row in rows]))
        else:
            self.assertTrue(True)

    def test02(self):
        if not self.xlrd_is_installed:
            rows = taxtastic.utils.get_new_nodes(os.path.join(datadir,'new_taxa.xls'))
            self.assertRaises(AttributeError, next, rows)
        else:
            self.assertTrue(True)

    def test03(self):
        rows = list(taxtastic.utils.get_new_nodes(os.path.join(datadir,'new_taxa.csv')))
        self.check_parent_id(rows)
        self.check_children(rows)

    def test04(self):
        rows = list(taxtastic.utils.get_new_nodes(os.path.join(datadir,'new_taxa_mac.csv')))
        self.check_parent_id(rows)

class StatsFileParsingMixIn(object):
    """
    Base class for stats file parsers.

    Classes using this mixing should define a property returning the

    Requires two files in datadir: example.txt and example.txt.json with expected results.
    """
    # Example. Will look for example.txt.json for expected
    test_file_name = 'example.txt'

    def setUp(self):
        self.test_path = config.data_path(self.test_file_name)
        self.expected_path = config.data_path(self.test_file_name + '.json')

    @property
    def expected(self):
        with open(self.expected_path) as fp:
            return json.load(fp)

    def test_parse(self):
        with open(self.test_path) as fp:
            result = self.parse_func(fp)
            # Hack: Pass through json so unicode matches
            actual = json.loads(json.dumps(result))
            self.assertEqual(self.expected, actual)

    @property
    def parse_func(self):
        # Example
        raise ValueError("Override")

class FastTreeStatsMixin(StatsFileParsingMixIn):
    @property
    def parse_func(self):
        return taxtastic.utils.parse_fasttree


class FastTreeMissingGTRTestCase(FastTreeStatsMixin, unittest.TestCase):
    test_file_name = 'fasttree-missing-rates.txt'
    def test_parse(self):
        with open(self.test_path) as fp:
            self.assertRaises(taxtastic.utils.InvalidLogError,
                    self.parse_func, fp)

class FastTreeDNATestCase(FastTreeStatsMixin, unittest.TestCase):
    test_file_name = 'fastree_dna_stats.txt'

class PhyMLStatsMixIn(StatsFileParsingMixIn):
    @property
    def parse_func(self):
        return taxtastic.utils.parse_phyml

class PhyMLAminoAcidTestCase(PhyMLStatsMixIn, unittest.TestCase):
    test_file_name = 'phyml_aa_stats.txt'

class PhyMLDNATestCase(PhyMLStatsMixIn, unittest.TestCase):
    test_file_name = 'phyml_dna_stats.txt'

class RAxMLStatsMixIn(StatsFileParsingMixIn):
    @property
    def parse_func(self):
        return taxtastic.utils.parse_raxml

class RAxMLAminoAcidTestCase(RAxMLStatsMixIn, unittest.TestCase):
    test_file_name = 'RAxML_info.aa'

class RAxMLDNATestCase(RAxMLStatsMixIn, unittest.TestCase):
    test_file_name = 'RAxML_info.testNuc'
