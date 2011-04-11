#!/usr/bin/env python

import sys
import os
import unittest
import logging
import shutil
import sqlite3

import config
import taxonomy
from taxonomy.errors import OperationalError, IntegrityError

log = logging

module_name = os.path.split(sys.argv[0])[1].rstrip('.py')
outputdir = os.path.abspath(config.outputdir)
datadir = os.path.abspath(config.datadir)

def newdir(path):
    """
    Create a new directory "path", deleting an existing directory if
    necessary.
    """

    shutil.rmtree(path, ignore_errors = True)
    os.makedirs(path)
    
class TestFetchData(unittest.TestCase):
    def setUp(self):
        self.funcname = '_'.join(self.id().split('.')[-2:])
        self.outdir = os.path.join(outputdir, self.funcname)
        newdir(self.outdir)
        _, self.zfilename = os.path.split(taxonomy.ncbi.ncbi_data_url)
        
    def test01(self):
        zfile = os.path.join(self.outdir, self.zfilename)        
        fout, downloaded = taxonomy.ncbi.fetch_data(dest_dir=self.outdir)

        # file is downloaded the first time
        self.assertTrue(downloaded)
        self.assertTrue(os.path.isfile(fout))
        self.assertTrue(zfile == fout)

        # ... but not the second time
        fout, downloaded = taxonomy.ncbi.fetch_data(dest_dir=self.outdir)
        self.assertFalse(downloaded)

        # ... unless clobber = True
        fout, downloaded = taxonomy.ncbi.fetch_data(dest_dir=self.outdir, clobber=True)
        self.assertTrue(downloaded)

        
class TestDbconnect(unittest.TestCase):

    def setUp(self):
        self.funcname = '_'.join(self.id().split('.')[-2:])
        self.dbname = os.path.join(outputdir, self.funcname + '.db')
        log.info(self.dbname)

    def test01(self):
        with taxonomy.ncbi.db_connect(self.dbname, clobber = True) as con:
            cur = con.cursor()
            cur.execute('select name from sqlite_master where type = "table"')
            tables = set(i for j in cur.fetchall() for i in j) # flattened
            self.assertTrue(set(['nodes','names','merged','source']).issubset(tables))

        # connection object is returned is schema already exists
        con = taxonomy.ncbi.db_connect(self.dbname, clobber = False)
        
class TestLoadData(unittest.TestCase):

    def setUp(self):
        self.funcname = '_'.join(self.id().split('.')[-2:])
        self.zfile = os.path.join(outputdir, 'taxdmp.zip')
        # reuse this after the first download
        _ , downloaded = taxonomy.ncbi.fetch_data(dest_dir = outputdir)
        self.dbname = os.path.join(outputdir, self.funcname + '.db')
        try:
            os.remove(self.dbname)
        except OSError:
            pass
        
    def test01(self):
        maxrows = 10
        with taxonomy.ncbi.db_connect(self.dbname, clobber = False) as con:
            taxonomy.ncbi.db_load(con, self.zfile, maxrows = maxrows)
            cur = con.cursor()
            cur.execute('select * from names')
            self.assertTrue(len(list(cur.fetchall())) == maxrows)
            
        with taxonomy.ncbi.db_connect(self.dbname, clobber = True) as con:
            taxonomy.ncbi.db_load(con, self.zfile, maxrows = 10)
            cur = con.cursor()
            cur.execute('select * from names')
            self.assertTrue(len(list(cur.fetchall())) == maxrows)

            # shouldn't be able to load data a second time, so number
            # of rows should not change
            taxonomy.ncbi.db_load(
                con = con,
                archive = self.zfile,
                maxrows = 10)
            cur.execute('select * from names')
            self.assertTrue(len(list(cur.fetchall())) == maxrows)


            
    # def test02(self):
    #     with taxonomy.ncbi.db_connect(self.dbname, clobber = True) as con:
    #         taxonomy.ncbi.db_load(con, self.zfile)




