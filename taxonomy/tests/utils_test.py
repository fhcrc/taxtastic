#!/usr/bin/env python

import sys
import os
import unittest
import logging
import itertools
import sqlite3
import shutil
import time
import pprint

import config
import Taxonomy

log = logging

outputdir = os.path.abspath(config.outputdir)
datadir = os.path.abspath(config.datadir)

class TestReadSpreadsheet(unittest.TestCase):

    def setUp(self):
        self.funcname = '_'.join(self.id().split('.')[-2:])

    def tearDown(self):
        pass

    def test01(self):
        headers, rows = Taxonomy.utils.read_spreadsheet(
            os.path.join(datadir,'new_taxa.xls'))
        check = lambda val: isinstance(val, float)
        self.assertTrue(all([check(row['parent_id']) for row in rows]))

    def test02(self):
        headers, rows = Taxonomy.utils.read_spreadsheet(
            os.path.join(datadir,'new_taxa.xls'),
            fmts={'tax_id':'%i','parent_id':'%i'}
            )
        check = lambda val: isinstance(val, str) and '.' not in val
        self.assertTrue(all([check(row['parent_id']) for row in rows]))


class TestGetNewNodes(unittest.TestCase):

    def setUp(self):
        self.funcname = '_'.join(self.id().split('.')[-2:])

    def tearDown(self):
        pass

    def test01(self):
        rows = Taxonomy.utils.get_new_nodes(os.path.join(datadir,'new_taxa.xls'))
        check = lambda val: isinstance(val, str) and '.' not in val
        self.assertTrue(all([check(row['parent_id']) for row in rows]))

class TestWriteConfig(unittest.TestCase):
    def setUp(self):
        self.funcname = '_'.join(self.id().split('.')[-2:])
        self.fname = os.path.join(outputdir, self.funcname) + '.ini'

    def tearDown(self):
        pass

    def test01(self):

        options = {('sec1','opt1'):'val1',('sec1','opt2'):'val2',('sec2','opt3'):'val3',('sec2','opt4'):'val4'}
        Taxonomy.package.write_config(
            fname = self.fname,
            optdict = options,
            sections = dict(sec1=['opt1','opt2'], sec2=['opt3','opt4'])
            )

    def test02(self):

        options = {('sec1','opt1'):'val1',('sec1','opt2'):'val2',('sec2','opt3'):'val3'}
        Taxonomy.package.write_config(
            fname = self.fname,
            optdict = options,
            sections = dict(sec1=['opt1','opt2'], sec2=['opt3','opt4'])
            )

