#!/usr/bin/env python

import sys
import os
from os import path
import unittest
import logging
import shutil
import sqlite3

import config
import taxtastic
import taxtastic.ncbi
from taxtastic.errors import OperationalError, IntegrityError
from taxtastic.utils import mkdir, rmdir

log = logging

outputdir = config.outputdir
datadir = config.datadir
ncbi_master_db = config.ncbi_master_db
ncbi_data = config.ncbi_data
    
# class TestFetchData(unittest.TestCase):
#     def setUp(self):
#         self.funcname = '_'.join(self.id().split('.')[-2:])
#         self.outdir = os.path.join(outputdir, self.funcname)
#         newdir(self.outdir)
#         _, self.zfilename = os.path.split(taxtastic.ncbi.ncbi_data_url)
        
#     def test01(self):
#         zfile = os.path.join(self.outdir, self.zfilename)        
#         fout, downloaded = taxtastic.ncbi.fetch_data(dest_dir=self.outdir)

#         # file is downloaded the first time
#         self.assertTrue(downloaded)
#         self.assertTrue(os.path.isfile(fout))
#         self.assertTrue(zfile == fout)

#         # ... but not the second time
#         fout, downloaded = taxtastic.ncbi.fetch_data(dest_dir=self.outdir)
#         self.assertFalse(downloaded)

#         # ... unless clobber = True
#         fout, downloaded = taxtastic.ncbi.fetch_data(dest_dir=self.outdir, clobber=True)
#         self.assertTrue(downloaded)
        
class TestDbconnect(unittest.TestCase):

    def setUp(self):
        self.funcname = '_'.join(self.id().split('.')[-2:])

    def test01(self):
        with taxtastic.ncbi.db_connect(ncbi_master_db, clobber = False) as con:
            cur = con.cursor()
            cur.execute('select name from sqlite_master where type = "table"')
            tables = set(i for j in cur.fetchall() for i in j) # flattened
            self.assertTrue(set(['nodes','names','merged','source']).issubset(tables))

class TestBase(unittest.TestCase):

    def setUp(self):
        self.funcname = '.'.join(self.id().split('.')[-1:])
        self.outdir = path.join(outputdir, self.funcname)
        self.dbname = os.path.join(self.outdir, 'taxonomy.db')
        mkdir(self.outdir, clobber = True)
    
class TestLoadData(TestBase):

    maxrows = 10
                    
    def test01(self):
        # we should be starting from scratch
        self.assertFalse(path.isfile(self.dbname))

        with taxtastic.ncbi.db_connect(self.dbname) as con:
            taxtastic.ncbi.db_load(con, ncbi_data, maxrows = self.maxrows)
            cur = con.cursor()
            cur.execute('select * from names')
            self.assertTrue(len(list(cur.fetchall())) == self.maxrows)

        # test clobber argument
        with taxtastic.ncbi.db_connect(self.dbname, clobber = True) as con:
            taxtastic.ncbi.db_load(con, ncbi_data, maxrows = self.maxrows)
            cur = con.cursor()
            cur.execute('select * from names')
            self.assertTrue(len(list(cur.fetchall())) == self.maxrows)

            # shouldn't be able to load data a second time, so number
            # of rows should not change
            taxtastic.ncbi.db_load(
                con = con,
                archive = ncbi_data,
                maxrows = self.maxrows)
            cur.execute('select * from names')
            self.assertTrue(len(list(cur.fetchall())) == self.maxrows)




