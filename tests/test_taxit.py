import json
import logging
from os import path
import shutil

from taxtastic import refpkg

from . import config
from .config import TestScriptBase

log = logging

TestScriptBase.executable = path.join(path.dirname(__file__), '..', 'taxit')
TestScriptBase.outputdir = config.outputdir
TestScriptBase.taxdb = config.ncbi_master_db
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

        contents_json = path.join(self.packagename,'CONTENTS.json')
        self.assertTrue(path.exists(contents_json))

        with open(contents_json) as f:
            contents = json.load(f)

        self.assertEqual(contents['metadata']['format_version'], refpkg.FORMAT_VERSION)

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
        """ Minimal test: a existing database is opened and that's it."""
        self.cmd_ok('taxtable -d %(taxdb)s > %(outfile)s')
        self.assertTrue(path.isfile(self.outfile))

    def test03(self):
        """ Minimal test: a existing database is opened and that's it."""
        self.cmd_ok('taxtable -d %(taxdb)s -o %(outfile)s')
        self.assertTrue(path.isfile(self.outfile))

    def test04(self):
        """Specify a single tax_id"""
        self.cmd_ok('taxtable -d %(taxdb)s -o %(outfile)s -t 180164')
        self.assertTrue(path.isfile(self.outfile))

    def test05(self):
        """Specify more than one tax_id"""
        self.cmd_ok('taxtable -d %(taxdb)s -o %(outfile)s -t 180164,166486')
        self.assertTrue(path.isfile(self.outfile))

    def test06(self):
        """taxids using an input file"""
        self.cmd_ok('taxtable -d %(taxdb)s -o %(outfile)s -t %(datadir)s/taxids1.txt')
        self.assertTrue(path.isfile(self.outfile))


