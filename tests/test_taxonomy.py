#!/usr/bin/env python

import os
from os import path
import logging
import shutil
import unittest

from sqlalchemy import create_engine

from . import config
from .config import TestBase

import taxtastic
from taxtastic.taxonomy import Taxonomy
import taxtastic.ncbi
import taxtastic.utils

log = logging

datadir = config.datadir

echo = False

zfile = config.ncbi_data
dbname = config.ncbi_master_db

try:
    import xlrd
except ImportError:
    xlrd = None

class TestTaxonomyBase(TestBase):

    def setUp(self):
        self.engine = create_engine('sqlite:///%s' % self.dbname, echo=echo)
        self.tax = Taxonomy(self.engine, taxtastic.ncbi.ranks)

    def tearDown(self):
        self.engine.dispose()


class TestAddNode(TestTaxonomyBase):
    """
    test tax.add_node
    """

    def setUp(self):
        self.dbname = path.join(self.mkoutdir(), 'taxonomy.db')
        log.info(self.dbname)
        shutil.copyfile(dbname, self.dbname)
        super(TestAddNode, self).setUp()

    def tearDown(self):
        pass

    def test01(self):
        self.tax.add_node(
            tax_id = '1578_1',
            parent_id = '1578',
            rank = 'species_group',
            tax_name = 'Lactobacillus helveticis/crispatus',
            source_id = 2
            )

        lineage = self.tax.lineage('1578_1')
        self.assertTrue(lineage['tax_id'] == '1578_1')
        self.assertTrue(lineage['tax_name'] == 'Lactobacillus helveticis/crispatus')

    def test02(self):

        new_taxid = '1578_1'
        new_taxname = 'Lactobacillus helveticis/crispatus'
        children = ['47770', # crispatus
                    '1587'] # helveticus

        self.tax.add_node(
            tax_id = new_taxid,
            parent_id = '1578',
            rank = 'species_group',
            tax_name = new_taxname,
            children = children,
            source_id = 2
            )

        lineage = self.tax.lineage(new_taxid)
        self.assertTrue(lineage['tax_id'] == new_taxid)
        self.assertTrue(lineage['tax_name'] == new_taxname)

        for taxid in children:
            lineage = self.tax.lineage(taxid)
            self.assertTrue(lineage['parent_id'] == new_taxid)

    @unittest.skipIf(xlrd is None, "xlrd is not installed")
    def test03(self):

        rows = taxtastic.utils.get_new_nodes(os.path.join(datadir,'new_taxa.xls'))
        for d in rows:
            d['source_id'] = 2
            self.tax.add_node(**d)

        new_taxid = '1578_1'
        new_taxname = 'Lactobacillus helveticis/crispatus'
        children = ['47770', # crispatus
                    '1587'] # helveticus
        lineage = self.tax.lineage(new_taxid)

        self.assertTrue(lineage['tax_id'] == new_taxid)
        self.assertTrue(lineage['tax_name'] == new_taxname)

        for taxid in children:
            lineage = self.tax.lineage(taxid)
            self.assertTrue(lineage['parent_id'] == new_taxid)
