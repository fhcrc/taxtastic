#!/usr/bin/env python

import re
import os
from os import path
import logging

import sqlalchemy as sa

import taxtastic
import taxtastic.ncbi
from taxtastic.ncbi import read_names, read_archive

from . import config
from .config import TestBase

log = logging

outputdir = config.outputdir
datadir = config.datadir
ncbi_master_db = config.ncbi_master_db
ncbi_data = config.ncbi_data


class TestDbconnect(TestBase):

    def test01(self):
        engine = sa.create_engine('sqlite:///' + ncbi_master_db)
        taxtastic.ncbi.db_connect(engine)
        with engine.begin() as con:
            result = con.execute(sa.text(
                'select name from sqlite_master where type = "table"'))
            tables = set(i[0] for i in result)
            self.assertTrue(
                set(['nodes', 'names', 'merged', 'source']).issubset(tables))


class TestLoadData(TestBase):

    def setUp(self):
        outdir = self.mkoutdir()
        self.db_path = os.path.join(outdir, 'taxonomy.db')
        self.url = 'sqlite:///' + self.db_path

    def test01(self):
        # we should be starting from scratch
        self.assertFalse(path.isfile(self.db_path))
        engine = sa.create_engine(self.url)

        taxtastic.ncbi.db_connect(engine)
        self.assertTrue(path.isfile(self.db_path))


class TestReadNames(TestBase):

    def setUp(self):
        self.zipfile = ncbi_data

    def test02(self):
        """
        is_classified always None
        """

        rows = read_names(rows=read_archive(self.zipfile, 'names.dmp'))
        headers = next(rows)
        is_classified = headers.index('is_classified')
        self.assertEqual(
            set(row[is_classified] for row in rows), set([None]))


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
