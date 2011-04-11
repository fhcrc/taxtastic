#!/usr/bin/env python


import os
import unittest
import logging

import config
import taxonomy

log = logging

outputdir = os.path.abspath(config.outputdir)
datadir = os.path.abspath(config.datadir)

if hasattr(taxonomy.utils, 'read_spreadsheet'):
    class TestReadSpreadsheet(unittest.TestCase):

        def setUp(self):
            self.funcname = '_'.join(self.id().split('.')[-2:])

        def tearDown(self):
            pass

        def test01(self):
            headers, rows = taxonomy.utils.read_spreadsheet(
                os.path.join(datadir,'new_taxa.xls'))
            check = lambda val: isinstance(val, float)
            self.assertTrue(all([check(row['parent_id']) for row in rows]))

        def test02(self):
            headers, rows = taxonomy.utils.read_spreadsheet(
                os.path.join(datadir,'new_taxa.xls'),
                fmts={'tax_id':'%i','parent_id':'%i'}
                )
            check = lambda val: isinstance(val, str) and '.' not in val
            self.assertTrue(all([check(row['parent_id']) for row in rows]))

## TODO: need to test creation of new nodes with csv (as opposed to xls) output
class TestGetNewNodes(unittest.TestCase):

    def setUp(self):
        self.funcname = '_'.join(self.id().split('.')[-2:])

    def tearDown(self):
        pass

    if hasattr(taxonomy.utils, 'read_spreadsheet'):
        def test01(self):
            rows = taxonomy.utils.get_new_nodes(os.path.join(datadir,'new_taxa.xls'))
            check = lambda val: isinstance(val, str) and '.' not in val
            self.assertTrue(all([check(row['parent_id']) for row in rows]))
    else:
        def test02(self):
            self.assertRaises(AttributeError, taxonomy.utils.get_new_nodes,
                              os.path.join(datadir,'new_taxa.xls'))

