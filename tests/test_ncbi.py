#!/usr/bin/env python

import os
from os import path
import logging
from itertools import groupby

import taxtastic
import taxtastic.ncbi
from taxtastic.ncbi import read_names, read_archive, UNCLASSIFIED_REGEX

from . import config
from .config import TestBase

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

class TestDbconnect(TestBase):

    def test01(self):
        with taxtastic.ncbi.db_connect(ncbi_master_db) as con:
            cur = con.cursor()
            cur.execute('select name from sqlite_master where type = "table"')
            tables = set(i for j in cur.fetchall() for i in j) # flattened
            self.assertTrue(set(['nodes','names','merged','source']).issubset(tables))

class TestLoadData(TestBase):

    maxrows = 10

    def setUp(self):
        outdir = self.mkoutdir()
        self.dbname = os.path.join(outdir, 'taxonomy.db')

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


class TestReadNames(TestBase):

    def setUp(self):
        self.zipfile = ncbi_data
        
    def test01(self):
        """
        is_classified always 0 or 1 if unclassified_regex is provided
        """
        
        rows = read_names(rows = read_archive(self.zipfile, 'names.dmp'),
                          unclassified_regex = UNCLASSIFIED_REGEX)
        self.assertEquals(set(row[-1] for row in rows), set([0,1]))

    def test02(self):
        """
        is_classified always None if unclassified_regex not provided
        """

        rows = read_names(rows = read_archive(self.zipfile, 'names.dmp'))
        self.assertEquals(set(row[-1] for row in rows), set([None]))


# class TestReadNamesExhaustively(TestReadNames):

#     def setUp(self):
#         self.zipfile, downloaded = taxtastic.ncbi.fetch_data(dest_dir = config.outputdir)

    
#     def test03(self):        
#         """
#         Print classified names
#         """
        
#         rows = read_names(rows = read_archive(self.zipfile, 'names.dmp'),
#                           unclassified_regex = UNCLASSIFIED_REGEX)

#         for row in rows:
#             if row[-2] and row[-1]:
#                 print row[1]

                
#     def test04(self):        
#         """
#         Print unclassified names
#         """
        
#         rows = read_names(rows = read_archive(self.zipfile, 'names.dmp'),
#                           unclassified_regex = UNCLASSIFIED_REGEX)

#         for row in rows:
#             if row[-2] and not row[-1]:
#                 print row[1]



