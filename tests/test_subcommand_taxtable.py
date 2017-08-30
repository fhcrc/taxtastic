#!/usr/bin/env python

from os import path
import logging

from sqlalchemy import create_engine

from . import config
from .config import TestBase

from taxtastic.taxonomy import Taxonomy

log = logging

datadir = config.datadir

echo = False

dbname = config.ncbi_master_db


class TaxTableSetup(TestBase):

    def setUp(self):
        self.engine = create_engine('sqlite:///%s' % dbname, echo=echo)
        self.tax = Taxonomy(self.engine)

    def tearDown(self):
        self.engine.dispose()


class TestTaxonomyInit(TaxTableSetup):

    def test01(self):
        self.tax._node('2')

    def test02(self):
        self.assertRaises(ValueError, self.tax._node, 'buh')


class TestGetLineagePrivate(TaxTableSetup):

    def test01(self):
        lineage = self.tax._get_lineage('1')
        self.assertTrue(lineage == [('root', '1')])

    def test02(self):
        tax_id = '1280'  # staph aureus

        lineage = self.tax._get_lineage(tax_id)
        self.assertTrue(lineage[0][0] == 'root')
        self.assertTrue(lineage[-1][0] == 'species')

    def test03(self):
        tax_id = '30630'  # deprecated; Microtus levis Taxonomy ID: 537919
        self.assertRaises(
            ValueError, self.tax._get_lineage, tax_id, merge_obsolete=False)

    def test04(self):
        tax_id = 'foo'
        self.assertRaises(ValueError, self.tax._get_lineage, tax_id)


class TestGetMerged(TaxTableSetup):

    def test01(self):
        tax_id = '1378'
        merged = self.tax._get_merged(tax_id)
        self.assertEqual(tax_id, merged)

    def test02(self):
        tax_id = '30630'  # deprecated; Microtus levis Taxonomy ID: 537919
        merged = self.tax._get_merged(tax_id)
        self.assertFalse(merged is None)


class TestTaxNameSearch(TaxTableSetup):

    def test01(self):
        tax_id, tax_name, is_primary = self.tax.primary_from_name('Gemella')
        self.assertTrue(tax_id == '1378')
        self.assertTrue(is_primary)

    def test02(self):
        self.assertRaises(ValueError, self.tax.primary_from_name, 'buggabugga')

    def test03(self):
        tax_id, tax_name, is_primary = self.tax.primary_from_name(
            'Gemella Berger 1960')

        self.assertTrue(tax_id == '1378')
        self.assertFalse(is_primary)


class TestSynonyms(TaxTableSetup):

    def test01(self):
        self.tax.synonyms(tax_id='1378')

    def test02(self):
        self.tax.synonyms(tax_name='Gemella')


class TestGetLineagePublic(TaxTableSetup):

    def test01(self):
        lineage = self.tax.lineage('1')

        self.assertTrue(lineage['root'] == '1')
        self.assertTrue(lineage['rank'] == 'root')

    def test02(self):
        tax_id = '1280'  # staph aureus
        lineage = self.tax.lineage(tax_id)
        self.assertTrue(lineage['rank'] == 'species')

        keys = set(lineage.keys())
        ranks = set(self.tax.ranks)
        self.assertTrue(
            keys - ranks == set(['parent_id', 'tax_id', 'rank', 'tax_name']))

    def test03(self):
        tax_id = '1378'  # Gemella; lineage has two successive no_rank taxa
        lineage = self.tax.lineage(tax_id)
        self.assertTrue(lineage['rank'] == 'genus')

        keys = set(lineage.keys())
        ranks = set(self.tax.ranks)
        self.assertTrue(
            keys - ranks == set(['parent_id', 'tax_id', 'rank', 'tax_name']))

    def test04(self):
        self.assertRaises(ValueError, self.tax.lineage,
                          tax_id=None, tax_name=None)

    def test05(self):
        self.assertRaises(ValueError, self.tax.lineage,
                          tax_id='1', tax_name='root')

    def test06(self):
        tax_id = '1378'  # Gemella; lineage has two successive no_rank taxa
        tax_name = 'Gemella'
        self.tax.lineage(tax_name=tax_name)

        lineage = self.tax.lineage(tax_id)
        self.assertEqual(lineage['rank'], 'genus')


class TestMethods(TaxTableSetup):

    def test01(self):
        taxname = self.tax.primary_from_id('1280')
        self.assertTrue(taxname == 'Staphylococcus aureus')

    def test02(self):
        self.assertRaises(ValueError, self.tax.primary_from_id, 'buh')

