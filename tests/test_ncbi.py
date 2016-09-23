#!/usr/bin/env python

import re
import os
from os import path
import logging

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


class TestDbconnect(TestBase):

    def test01(self):
        engine = taxtastic.ncbi.db_connect(ncbi_master_db)
        with engine.begin() as con:
            result = con.execute(
                'select name from sqlite_master where type = "table"')
            tables = set(i[0] for i in result)
            self.assertTrue(
                set(['nodes', 'names', 'merged', 'source']).issubset(tables))


class TestLoadData(TestBase):

    names_rows_count = 1100

    def setUp(self):
        outdir = self.mkoutdir()
        self.dbname = os.path.join(outdir, 'taxonomy.db')

    def test01(self):
        # we should be starting from scratch
        self.assertFalse(path.isfile(self.dbname))

        engine = taxtastic.ncbi.db_connect(self.dbname)
        taxtastic.ncbi.db_load(engine, ncbi_data)
        with engine.begin() as conn:
            result = conn.execute('select 1 AS i from names')
            self.assertEqual(self.names_rows_count, len(list(result)))

        # test clobber argument
        engine = taxtastic.ncbi.db_connect(self.dbname, clobber=True)
        taxtastic.ncbi.db_load(engine, ncbi_data)
        with engine.begin() as conn:
            result = conn.execute('select 1 AS i from names')
            self.assertEqual(self.names_rows_count, len(list(result)))

            # shouldn't be able to load data a second time
            with self.assertRaises(taxtastic.errors.IntegrityError):
                taxtastic.ncbi.db_load(engine, archive=ncbi_data)


class TestReadNames(TestBase):

    def setUp(self):
        self.zipfile = ncbi_data

    def test01(self):
        """
        is_classified always 0 or 1 if unclassified_regex is provided
        """

        rows = read_names(rows=read_archive(self.zipfile, 'names.dmp'),
                          unclassified_regex=UNCLASSIFIED_REGEX)
        self.assertEquals(set(row['is_classified']
                              for row in rows), set([0, 1]))

    def test02(self):
        """
        is_classified always None if unclassified_regex not provided
        """

        rows = read_names(rows=read_archive(self.zipfile, 'names.dmp'))
        self.assertEquals(set(row['is_classified']
                              for row in rows), set([None]))


class TestUnclassifiedRegex(TestBase):
    """
    Test the heuristic used to determine if a taxonomic name is meaningful.
    """

    def setUp(self):
        self.pieces = taxtastic.ncbi.UNCLASSIFIED_REGEX_COMPONENTS
        self.regexes = [re.compile(piece) for piece in self.pieces]
        with open(config.data_path('type_strain_names.txt')) as fp:
            self.type_strain_names = [i.rstrip() for i in fp]

    def test_no_type_strains_match(self):
        for strain_name in self.type_strain_names:
            for regex in self.regexes:
                m = regex.search(strain_name)
                if m:
                    self.fail('"{0}" matches "{1}"'.format(
                        strain_name, regex.pattern))

# def generate_test_unclassified_regex():
    #"""
    # Generate a test class verifying that none of the type strains in
    # type_strain_names.txt match the unclassified regex.
    #"""
    # def generate_test(strain_name):
        # def do_test(self):
            # for regex in self.regexes:
                #m = regex.search(strain_name)
                # if m:
                    #self.fail('"{0}" matches "{1}"'.format(strain_name, regex.pattern))
        # return do_test

    # class TestUnclassifiedRegex(TestBase):
        # def setUp(self):
            #self.pieces = taxtastic.ncbi.UNCLASSIFIED_REGEX_COMPONENTS
            #self.regexes = [re.compile(piece) for piece in self.pieces]

    # with open(config.data_path('type_strain_names.txt')) as fp:
        #type_strain_names = [i.rstrip() for i in fp]

    # for s in type_strain_names:
        #test_fn = generate_test(s)
        #func_name = 'test_{0}_no_match'.format(re.sub(r'[ -.]', '_', s).lower())
        #test_fn.__name__ = func_name
        #setattr(TestUnclassifiedRegex, func_name, test_fn)

    # return TestUnclassifiedRegex

#TestUnclassifiedRegex = generate_test_unclassified_regex()
#del generate_test_unclassified_regex
