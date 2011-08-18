#!/usr/bin/env python


import os
import unittest
import logging
import pprint

import config
import taxtastic
import taxtastic.utils


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
            self.assertRaises(
                AttributeError, taxtastic.utils.get_new_nodes,
                os.path.join(datadir,'new_taxa.xls'))
        else:
            self.assertTrue(True)

    def test03(self):
        rows = list(taxtastic.utils.get_new_nodes(os.path.join(datadir,'new_taxa.csv')))
        self.check_parent_id(rows)
        self.check_children(rows)

    def test04(self):
        rows = list(taxtastic.utils.get_new_nodes(os.path.join(datadir,'new_taxa_mac.csv')))
        self.check_parent_id(rows)
        


