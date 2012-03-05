#!/usr/bin/env python

import os
from os import path
import logging
import shutil
import unittest

from sqlalchemy import create_engine

import config
from config import TestBase

import taxtastic
from taxtastic.taxonomy import Taxonomy, ranks_below, is_below
import taxtastic.ncbi
import taxtastic.utils

log = logging

datadir = config.datadir

echo = False

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

def test__node():
    engine = create_engine('sqlite:///../testfiles/small_taxonomy.db', echo=False)
    tax = Taxonomy(engine, taxtastic.ncbi.ranks)
    assert tax._node(None) == None
    assert tax._node('91061') == (u'1239', u'class')

def test_sibling_of():
    engine = create_engine('sqlite:///../testfiles/taxonomy.db', echo=False)
    tax = Taxonomy(engine, taxtastic.ncbi.ranks)
    assert tax.sibling_of(None) == None
    assert tax.sibling_of('91061') == '186801'
    assert tax.sibling_of('1696') == None

def test_child_of():
    engine = create_engine('sqlite:///../testfiles/small_taxonomy.db', echo=False)
    tax = Taxonomy(engine, taxtastic.ncbi.ranks)
    assert tax.child_of(None) == None
    assert tax.child_of('1239') == '91061'
    assert tax.children_of('1239', 2) == ['91061', '186801']
    
def test_is_ancestor_of():
    engine = create_engine('sqlite:///../testfiles/taxonomy.db', echo=False)
    tax = Taxonomy(engine, taxtastic.ncbi.ranks)
    assert tax.is_ancestor_of('1280','1239')
    assert tax.is_ancestor_of(None, '1239') == False
    assert tax.is_ancestor_of('1239', None) == False

def test_rank_and_parent():
    engine = create_engine('sqlite:///../testfiles/taxonomy.db', echo=False)
    tax = Taxonomy(engine, taxtastic.ncbi.ranks)
    assert tax.rank(None) == None
    assert tax.rank('1239') == 'phylum'
    assert tax.rank('1280') == 'species'
    assert tax.parent_id(None) == None
    assert tax.parent_id('1239') == '2'

def test_species_below():
    engine = create_engine('sqlite:///../testfiles/taxonomy.db', echo=False)
    tax = Taxonomy(engine, taxtastic.ncbi.ranks)
    t = tax.species_below('1239')
    parent_id, rank = tax._node(t)
    for t in [None, '1239','186801','1117']:
        s = tax.species_below(t)
        assert t is None or s is None or tax.is_ancestor_of(s, t)
        assert s is None or tax.rank(s) == 'species'

def test_is_below():
    assert is_below('species', 'family')
    assert is_below('family', 'kingdom')
    assert not(is_below('kingdom', 'family'))
    assert ranks_below('species') == []
    assert ranks_below('family') == ['species','genus']

def test_nary_subtree():
    engine = create_engine('sqlite:///../testfiles/small_taxonomy.db', echo=False)
    tax = Taxonomy(engine, taxtastic.ncbi.ranks)
    assert tax.nary_subtree(None) == None
    t = tax.nary_subtree('1239')
    assert t == ['1280', '372074', '1579', '1580', '37734', '420335', '166485', '166486']
    
