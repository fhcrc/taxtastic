import json
import logging
from os import path
import shutil
import sys

from taxtastic.scripts.taxit import main
from taxtastic import refpkg

from . import config
from .config import TestScriptBase, TestBase

log = logging

TestScriptBase.executable = path.join(path.dirname(__file__), '..', 'taxit.py')
TestScriptBase.outputdir = config.outputdir
TestScriptBase.taxdb = 'sqlite:///' + config.ncbi_master_db
TestScriptBase.datadir = config.datadir


class TestHelp(TestScriptBase):

    def test01(self):
        self.cmd_ok('-h')

    def test02(self):
        self.cmd_ok('--help')

    def test03(self):
        self.cmd_fails('')

    def test04(self):
        self.cmd_fails('notacommand')

    def test05(self):
        self.cmd_ok('help create')

    def test06(self):
        self.cmd_ok('help check')

    def test07(self):
        self.cmd_ok('help help')

    def test08(self):
        self.cmd_fails('help notacommand')

    def test09(self):
        self.cmd_ok('--version')


class TestCreate(TestScriptBase):

    def setUp(self):
        super(TestCreate, self).setUp()
        self.packagename = path.join(self.mkoutdir(), 'refpkg')
        if path.exists(self.packagename):
            shutil.rmtree(self.packagename)

    def test01(self):
        self.cmd_fails('create -P %(packagename)s')

    def test02(self):
        """
        Create a minimal package.
        """

        self.cmd_ok('create -P %(packagename)s -l 16s')
        self.assertTrue(path.exists(self.packagename))

        contents_json = path.join(self.packagename, 'CONTENTS.json')
        self.assertTrue(path.exists(contents_json))

        with open(contents_json) as f:
            contents = json.load(f)

        self.assertEqual(contents['metadata'][
                         'format_version'], refpkg.FORMAT_VERSION)

        # test the --clobber option
        self.cmd_ok('create -P %(packagename)s -l 16s --clobber')

        # fails without --clobber because package already exists
        self.cmd_fails('create -P %(packagename)s -l 16s')


class TestTaxTable(TestScriptBase):
    """
    Unit tests for the taxtable sub-command.
    """

    def setUp(self):
        super(TestTaxTable, self).setUp()
        self.outfile = path.join(self.mkoutdir(), 'taxtable.csv')

    def test01(self):
        """
        Invalid arguments should cause a failure.
        """
        self.cmd_fails('taxtable --not-an-argument')
        self.assertFalse(path.isfile(self.outfile))

    def test02(self):
        """Specify a single tax_id"""
        self.cmd_ok('taxtable %(taxdb)s -o %(outfile)s -t 1280')
        self.assertTrue(path.isfile(self.outfile))

    def test03(self):
        """Specify more than one tax_id"""
        self.cmd_ok('taxtable %(taxdb)s -o %(outfile)s -t 1280 1281')
        self.assertTrue(path.isfile(self.outfile))

    def test04(self):
        """taxids using an input file"""
        self.cmd_ok(
            'taxtable %(taxdb)s -o %(outfile)s -f %(datadir)s/taxids1.txt')
        self.assertTrue(path.isfile(self.outfile))


class LonelyNodesTestCase(TestScriptBase):

    def setUp(self):
        super(LonelyNodesTestCase, self).setUp()
        self.outfile = path.join(self.mkoutdir(), 'lonely.txt')
        self.refpkg = config.data_path('lactobacillus2-0.2.refpkg')

    def test_all_ranks(self):
        self.cmd_ok('lonelynodes %(refpkg)s -o %(outfile)s')
        self.assertTrue(path.isfile(self.outfile))
        with open(self.outfile, **self.openargs) as fp:
            self.assertEqual("""tax_name,tax_id,rank
Bacilli,91061,class
Bacteria,2,domain
Enterobacteriaceae,543,family
Enterobacteriales,91347,order
Escherichia,561,genus
Escherichia coli,562,species
Gammaproteobacteria,1236,class
Lactobacillaceae,33958,family
Lactobacillales,186826,order
Lactobacillus,1578,genus
cellular organisms,131567,cellular_root
""", fp.read())

    def test_species(self):
        self.cmd_ok('lonelynodes %(refpkg)s -o %(outfile)s -r species')
        self.assertTrue(path.isfile(self.outfile))
        with open(self.outfile, **self.openargs) as fp:
            line = 'tax_name,tax_id,rank\nEscherichia coli,562,species'
            self.assertEqual(line, fp.read().strip())
