#!/usr/bin/env python

import os
from os import path
import logging
import shutil

from sqlalchemy import create_engine

from . import config
from .config import TestBase

import taxtastic
from taxtastic.taxonomy import Taxonomy, TaxonIntegrityError
import taxtastic.ncbi
import taxtastic.utils

log = logging

datadir = config.datadir

echo = False

dbname = config.ncbi_master_db


class TestTaxonomyBase(TestBase):

    def setUp(self):
        self.engine = create_engine('sqlite:///' + self.dbname, echo=echo)
        self.tax = Taxonomy(self.engine, taxtastic.ncbi.RANKS)

    def tearDown(self):
        self.engine.dispose()


class TestAddNode(TestTaxonomyBase):

    def setUp(self):
        self.dbname = path.join(self.mkoutdir(), 'taxonomy.db')
        log.info(self.dbname)
        shutil.copyfile(dbname, self.dbname)
        super(TestAddNode, self).setUp()

    def tearDown(self):
        pass

    def test01(self):
        self.tax.add_node(
            tax_id='1280_1',
            parent_id='1280',
            rank='subspecies',
            names=[{'tax_name': 'foo'}],
            source_name='ncbi'
        )

        lineage = self.tax.lineage('1280_1')
        self.assertEqual(lineage['tax_id'], '1280_1')
        self.assertEqual(lineage['tax_name'], 'foo')

    def test02(self):

        new_taxid = '1279_1'
        new_taxname = 'between genus and species'
        children = ['1280', '1281']

        self.tax.add_node(
            tax_id=new_taxid,
            parent_id='1279',
            rank='species_group',
            names=[{'tax_name': new_taxname}],
            children=children,
            source_name='foo'
        )

        lineage = self.tax.lineage(new_taxid)
        self.assertTrue(lineage['tax_id'] == new_taxid)
        self.assertTrue(lineage['tax_name'] == new_taxname)

        for taxid in children:
            lineage = self.tax.lineage(taxid)
            self.assertTrue(lineage['parent_id'] == new_taxid)

    def test03(self):

        new_taxid = '1279_1'
        new_taxname = 'between genus and species'
        children = ['1280', '1281']

        self.assertRaises(
            TaxonIntegrityError,
            self.tax.add_node,
            tax_id=new_taxid,
            parent_id='1279',
            rank='genus',
            names=[{'tax_name': new_taxname}],
            children=children,
            source_name='ncbi')

    def test04(self):
        # existing node
        self.assertRaises(
            ValueError,
            self.tax.add_node,
            tax_id='1280',
            parent_id='1279',
            rank='species',
            names=[{'tax_name': 'I already exist'}],
            source_name='ncbi'
        )

    def test05(self):
        self.tax.add_node(
            tax_id='1280_1',
            parent_id='1280',
            rank='subspecies',
            names=[
                {'tax_name': 'foo', 'is_primary': True},
                {'tax_name': 'bar'},
            ],
            source_name='ncbi'
        )

        lineage = self.tax.lineage('1280_1')
        self.assertEqual(lineage['tax_id'], '1280_1')
        self.assertEqual(lineage['tax_name'], 'foo')

    def test06(self):
        # multiple names, none primary
        self.assertRaises(
            ValueError,
            self.tax.add_node,
            tax_id='1280_1',
            parent_id='1280',
            rank='subspecies',
            names=[
                {'tax_name': 'foo'},
                {'tax_name': 'bar'},
            ],
            source_name='ncbi')

    def test07(self):
        self.tax.add_node(
            tax_id='1280_1',
            parent_id='1280',
            rank='subspecies',
            names=[
                {'tax_name': 'foo', 'is_primary': True},
                {'tax_name': 'bar'},
            ],
            source_name='ncbi',
            execute=False
        )

        self.assertRaises(ValueError, self.tax.lineage, '1280_1')

    def test08(self):
        # test has_node()
        self.assertTrue(self.tax.has_node('1280'))
        self.assertFalse(self.tax.has_node('foo'))


class TestAddName(TestTaxonomyBase):
    """
    test tax.add_node
    """

    def count_names(self, tax_id):
        with self.tax.engine.connect() as con:
            result = con.execute(
                'select count(*) from names where tax_id = ?', (tax_id,))
            return result.fetchone()[0]

    def count_primary_names(self, tax_id):
        with self.tax.engine.connect() as con:
            result = con.execute(
                'select count(*) from names where tax_id = ? and is_primary',
                (tax_id,))
            return result.fetchone()[0]

    def primary_name(self, tax_id):
        with self.tax.engine.connect() as con:
            result = con.execute(
                'select tax_name from names where tax_id = ? and is_primary',
                (tax_id,))
            val = result.fetchone()
            return val[0] if val else None

    def setUp(self):
        self.dbname = path.join(self.mkoutdir(), 'taxonomy.db')
        log.info(self.dbname)
        shutil.copyfile(dbname, self.dbname)
        super(TestAddName, self).setUp()

    def test_name01(self):
        names_before = self.count_names('1280')
        self.tax.add_name(tax_id='1280', tax_name='SA', source_name='ncbi')
        self.assertEqual(names_before + 1, self.count_names('1280'))

    def test_name02(self):
        # number of primary names should remain 1
        names_before = self.count_names('1280')
        self.assertEqual(self.count_primary_names('1280'), 1)
        self.tax.add_name(tax_id='1280', tax_name='SA', is_primary=True,
                          source_name='ncbi')
        self.tax.add_name(tax_id='1280', tax_name='SA2', is_primary=True,
                          source_name='ncbi')
        self.assertEqual(names_before + 2, self.count_names('1280'))
        self.assertEqual(self.count_primary_names('1280'), 1)

    def test_name03(self):
        # insertion of duplicate row fails
        self.tax.add_name(tax_id='1280', tax_name='SA', is_primary=True,
                          source_name='ncbi')

        self.assertRaises(
            ValueError, self.tax.add_name, tax_id='1280', tax_name='SA',
            is_primary=True, source_name='ncbi')

        self.assertEqual(self.primary_name('1280'), 'SA')


class TestGetSource(TestTaxonomyBase):
    def setUp(self):
        self.dbname = dbname
        super(TestGetSource, self).setUp()

    def test01(self):
        self.assertRaises(ValueError, self.tax.get_source)

    def test02(self):
        self.assertRaises(ValueError, self.tax.get_source, 1, 'ncbi')

    def test03(self):
        result = self.tax.get_source(source_id=1)
        self.assertDictEqual(result, {
            'description': 'ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdmp.zip',
            'id': 1, 'name': 'ncbi'})

    def test04(self):
        result = self.tax.get_source(source_name='ncbi')
        self.assertDictEqual(result, {
            'description': 'ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdmp.zip',
            'id': 1, 'name': 'ncbi'})

    def test05(self):
        self.assertRaises(ValueError, self.tax.get_source, source_id=2)


class TestAddSource(TestTaxonomyBase):
    def setUp(self):
        self.dbname = path.join(self.mkoutdir(), 'taxonomy.db')
        log.info(self.dbname)
        shutil.copyfile(dbname, self.dbname)
        super(TestAddSource, self).setUp()

    def tearDown(self):
        pass

    def sources(self):
        with self.tax.engine.connect() as con:
            result = con.execute('select * from source')
            return result.fetchall()

    def test01(self):
        self.tax.add_source('foo')
        self.assertEqual(self.sources()[1], (2, 'foo', None))

    def test02(self):
        self.tax.add_source('ncbi')
        self.assertEqual(
            self.sources(),
            [(1, 'ncbi', 'ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdmp.zip')])


def test__node():
    engine = create_engine(
        'sqlite:///../testfiles/small_taxonomy.db', echo=False)
    tax = Taxonomy(engine, taxtastic.ncbi.RANKS)
    assert tax._node(None) is None
    assert tax._node('91061') == ('1239', 'class')


def test_sibling_of():
    engine = create_engine('sqlite:///../testfiles/taxonomy.db', echo=False)
    tax = Taxonomy(engine, taxtastic.ncbi.RANKS)
    assert tax.sibling_of(None) is None
    assert tax.sibling_of('91061') == '186801'
    assert tax.sibling_of('1696') is None


def test_child_of():
    engine = create_engine(
        'sqlite:///../testfiles/small_taxonomy.db', echo=False)
    tax = Taxonomy(engine, taxtastic.ncbi.RANKS)
    assert tax.child_of(None) is None
    assert tax.child_of('1239') == '91061'
    assert tax.children_of('1239', 2) == ['91061', '186801']


def test_is_ancestor_of():
    engine = create_engine('sqlite:///../testfiles/taxonomy.db', echo=False)
    tax = Taxonomy(engine, taxtastic.ncbi.RANKS)
    assert tax.is_ancestor_of('1280', '1239')
    assert tax.is_ancestor_of(None, '1239') is False
    assert tax.is_ancestor_of('1239', None) is False


def test_rank_and_parent():
    engine = create_engine('sqlite:///../testfiles/taxonomy.db', echo=False)
    tax = Taxonomy(engine, taxtastic.ncbi.RANKS)
    assert tax.rank(None) is None
    assert tax.rank('1239') == 'phylum'
    assert tax.rank('1280') == 'species'
    assert tax.parent_id(None) is None
    assert tax.parent_id('1239') == '2'


def test_species_below():
    engine = create_engine('sqlite:///../testfiles/taxonomy.db', echo=False)
    tax = Taxonomy(engine, taxtastic.ncbi.RANKS)
    t = tax.species_below('1239')
    parent_id, rank = tax._node(t)
    for t in [None, '1239', '186801', '1117']:
        s = tax.species_below(t)
        assert t is None or s is None or tax.is_ancestor_of(s, t)
        assert s is None or tax.rank(s) == 'species'


def test_is_below():
    assert Taxonomy.is_below('species', 'family')
    assert Taxonomy.is_below('family', 'kingdom')
    assert not Taxonomy.is_below('kingdom', 'family')
    assert Taxonomy.ranks_below('species') == []
    assert Taxonomy.ranks_below('family') == ['species', 'genus']


def test_nary_subtree():
    engine = create_engine(
        'sqlite:///../testfiles/small_taxonomy.db', echo=False)
    tax = Taxonomy(engine, taxtastic.ncbi.RANKS)
    assert tax.nary_subtree(None) is None
    t = tax.nary_subtree('1239')
    assert t == ['1280', '372074', '1579', '1580',
                 '37734', '420335', '166485', '166486']
