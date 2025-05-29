#!/usr/bin/env python

from os import path
import logging
import shutil
from itertools import groupby

import sqlalchemy as sa
from sqlalchemy import create_engine

from . import config
from .config import TestBase, data_path

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
        self.tax = Taxonomy(self.engine)

    def tearDown(self):
        self.engine.dispose()


class TestAddNode(TestTaxonomyBase):

    def setUp(self):
        self.dbname = path.join(self.mkoutdir(), 'taxonomy.db')
        # log.info(self.dbname)
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
    test tax.add_name
    """

    def query(self, query, **params):
        with self.tax.engine.connect() as con:
            result = con.execute(sa.text(query), params)
            return result.fetchone()

    def count_names(self, tax_id):
        result = self.query(
            'select count(*) from names where tax_id=:tax_id', tax_id=tax_id)
        return result[0]

    def count_primary_names(self, tax_id):
        result = self.query(
            'select count(*) from names where tax_id=:tax_id and is_primary',
            tax_id=tax_id)
        return result[0]

    def primary_name(self, tax_id):
        result = self.query(
            'select tax_name from names where tax_id=:tax_id and is_primary',
            tax_id=tax_id)
        return result[0] if result else None

    def setUp(self):
        self.dbname = path.join(self.mkoutdir(), 'taxonomy.db')
        # log.info(self.dbname)
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


class TestAddNames(TestTaxonomyBase):
    """
    test tax.add_names
    """

    def setUp(self):
        self.dbname = path.join(self.mkoutdir(), 'taxonomy.db')
        # log.info(self.dbname)
        shutil.copyfile(dbname, self.dbname)
        super(TestAddNames, self).setUp()

    def test01(self):

        names = [
            dict(tax_id='1280', tax_name='SA', is_primary=True,
                 source_name='ncbi'),
            dict(tax_id='1280', tax_name='SA2', is_primary=False,
                 source_name='ncbi'),
        ]

        count_before = self.tax.fetchone(
            sa.text("select count(*) from names where tax_id = :tax_id"),
            tax_id='1280')[0]

        self.tax.add_names(tax_id='1280', names=names)

        count_after = self.tax.fetchone(
            sa.text("select count(*) from names where tax_id = :tax_id"),
            tax_id='1280')[0]

        self.assertEqual(count_before + 2, count_after)

    def test02(self):

        names = [
            dict(tax_id='1280', tax_name='SA', is_primary=True,
                 source_name='ncbi'),
            dict(tax_id='1280', tax_name='SA2', is_primary=True,
                 source_name='ncbi'),
        ]

        self.assertRaises(
            ValueError, self.tax.add_names, tax_id='1280', names=names)

    def test03(self):

        names = [
            dict(tax_id='1280', tax_name='SA', is_primary=True,
                 source_name='ncbi'),
            dict(tax_id='1280', tax_name='SA', is_primary=False,
                 source_name='ncbi'),
        ]

        self.assertRaises(
            ValueError, self.tax.add_names, tax_id='1280', names=names)


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
        # log.info(self.dbname)
        shutil.copyfile(dbname, self.dbname)
        super(TestAddSource, self).setUp()

    def tearDown(self):
        pass

    def sources(self):
        with self.tax.engine.connect() as con:
            result = con.execute(sa.text('select * from source'))
            return result.fetchall()

    def test01(self):
        res = self.tax.add_source('foo')
        self.assertEqual(res, (2, True))
        self.assertEqual(self.sources()[1], (2, 'foo', None))

    def test02(self):
        res = self.tax.add_source('ncbi')
        self.assertEqual(res, (1, False))
        self.assertEqual(
            self.sources(),
            [(1, 'ncbi', 'ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdmp.zip')])


class TestTaxonomyTree(TestTaxonomyBase):
    def setUp(self):
        self.dbname = data_path('small_taxonomy.db')
        super(TestTaxonomyTree, self).setUp()

    def tearDown(self):
        pass

    def test__node(self):
        self.assertRaises(ValueError, self.tax._node, None)
        self.assertEqual(self.tax._node('91061'), ('1239', 'class'))

    def test_sibling_of(self):
        self.assertRaises(ValueError, self.tax.sibling_of, None)
        self.assertEqual(self.tax.sibling_of('45669'), '1279')
        self.assertIsNone(self.tax.sibling_of('91061'))

    def test_child_of(self):
        self.assertRaises(ValueError, self.tax.child_of, None)
        self.assertEqual(self.tax.child_of('1239'), '91061')

    def test_children_of(self):
        self.assertEqual(self.tax.children_of('90964', 2), ['1279', '45669'])

    def test_is_ancestor_of(self):
        self.assertTrue(self.tax.is_ancestor_of('1280', '1239'))
        self.assertFalse(self.tax.is_ancestor_of(None, '1239'))
        self.assertFalse(self.tax.is_ancestor_of('1239', None))

    def test_rank_and_parent(self):
        self.assertRaises(ValueError, self.tax.rank, None)
        self.assertEqual(self.tax.rank('1239'), 'phylum')
        self.assertEqual(self.tax.rank('1280'), 'species')
        self.assertRaises(ValueError, self.tax.parent_id, None)
        self.assertEqual(self.tax.parent_id('1239'), '1783272')

    def test_species_below(self):
        t = self.tax.species_below('1239')
        parent_id, rank = self.tax._node(t)
        for t in [None, '1239', '186801', '1117']:
            s = self.tax.species_below(t)
            self.assertTrue(
                t is None or s is None or self.tax.is_ancestor_of(s, t))
            self.assertTrue(s is None or self.tax.rank(s) == 'species')

    def test_is_below(self):
        self.assertTrue(self.tax.is_below('species', 'family'))
        self.assertTrue(self.tax.is_below('family', 'kingdom'))
        self.assertFalse(self.tax.is_below('kingdom', 'family'))
        self.assertEqual(
            self.tax.ranks_below('species'),
            ['forma', 'subvariety', 'varietas', 'serogroup',
             'pathogroup', 'morph', 'genotype', 'biotype', 'subspecies'])
        self.assertEqual(
            self.tax.ranks_below('family'),
            ['forma', 'subvariety', 'varietas', 'serogroup', 'pathogroup',
             'morph', 'genotype', 'biotype', 'subspecies', 'species',
             'species_subgroup', 'species_group', 'subseries', 'series',
             'subsection', 'section', 'subgenus', 'genus', 'subtribe',
             'tribe', 'subfamily'])

    def test_nary_subtree(self):
        self.assertRaises(ValueError, self.tax.nary_subtree, None)
        self.assertEqual(
            self.tax.nary_subtree('1239'),
            ['1280', '1281', '45670', '138846'])


class TestGetLineageTable(TestTaxonomyBase):

    def setUp(self):
        self.dbname = path.join(self.mkoutdir(), 'taxonomy.db')
        # log.info(self.dbname)
        shutil.copyfile(dbname, self.dbname)
        super(TestGetLineageTable, self).setUp()

    def tearDown(self):
        pass

    def test01(self):
        lineage = self.tax._get_lineage_table(['1280', '1379'])

        species = {}
        for tax_id, grp in groupby(lineage, lambda row: getattr(row, 'tid')):
            species[tax_id] = list(grp)[-1]

        self.assertEqual(species['1280'].tax_name, 'Staphylococcus aureus')
        self.assertEqual(species['1379'].tax_name, 'Gemella haemolysans')


