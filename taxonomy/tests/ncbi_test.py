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

module_name = os.path.split(sys.argv[0])[1].rstrip('.py')
outputdir = os.path.abspath(config.outputdir)
datadir = os.path.abspath(config.datadir)

class TestDownload(unittest.TestCase):
    def setUp(self):
        self.funcname = '_'.join(self.id().split('.')[-2:])

    def test01(self):
        zfile = os.path.join(outputdir, 'taxdmp.zip')
        # os.remove(zfile)

        archive = Taxonomy.ncbi.fetch_data(dest_dir=outputdir)
        self.assertTrue(os.path.isfile(archive))
        self.assertTrue(zfile == archive)

        archive = Taxonomy.ncbi.fetch_data(dest_dir=outputdir)


class TestCreateSchema(unittest.TestCase):

    def setUp(self):
        self.funcname = '_'.join(self.id().split('.')[-2:])
        self.dbname = os.path.join(outputdir, self.funcname + '.db')
        log.info(self.dbname)

    def test01(self):
        con = Taxonomy.ncbi.db_connect(self.dbname, new=True)

class TestLoadData(unittest.TestCase):

    def setUp(self):
        self.funcname = '_'.join(self.id().split('.')[-2:])
        self.dbname = os.path.join(outputdir, self.funcname + '.db')
        self.zfile = os.path.join(outputdir, 'taxdmp.zip')

    def test01(self):
        con = Taxonomy.ncbi.db_connect(self.dbname, new=True)
        Taxonomy.ncbi.db_load(con, self.zfile, maxrows=10)
        con.close()

    def test02(self):
        con = Taxonomy.ncbi.db_connect(self.dbname, new=True)
        Taxonomy.ncbi.db_load(con, self.zfile, maxrows=10)
        con.close()



