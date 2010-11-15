#!/usr/bin/env python

import sys
import os
import unittest
import logging
import itertools
import sqlite3
import shutil
import time
import pprint
import commands

import config

log = logging

module_name = os.path.split(sys.argv[0])[1].rstrip('.py')
outputdir = os.path.abspath(config.outputdir)
datadir = os.path.abspath(config.datadir)

class TestBase(unittest.TestCase):

    executable = 'taxomatic.py'

    def __getitem__(self, i):
        """
        Enables string formatting, eg:
        print 'current function: %(funcname)s' % self
        """
        return getattr(self, i)

    def wrap_cmd(self, cmd=None, args=None):
        if cmd is None:
            cmd = self.executable
        input = cmd + ' ' + args
        log.info('\n >> '+ input)
        status, output = commands.getstatusoutput(input)
        log.info(output)
        return status, output

    def cmd_ok(self, cmd=None, args=None):
        status, output = self.wrap_cmd(cmd, args)
        self.assertTrue(status == 0)

    def cmd_fails(self, cmd=None, args=None):
        status, output = self.wrap_cmd(cmd, args)
        self.assertFalse(status == 0)

    def setUp(self):
        self.funcname = '_'.join(self.id().split('.')[-2:])
        loglevel = log.getLogger().level
        self.verbosity = {
            logging.WARNING:'',
            logging.INFO:'-v',
            logging.DEBUG:'-vv'}.get(loglevel,'')

        self.outputdir = config.outputdir
        self.outfile = os.path.join(self.outputdir, self.funcname)

    def tearDown(self):
        pass

class TestHelp(TestBase):

    def test01(self):
        self.cmd_ok(args='-h')

    def test02(self):
        self.cmd_ok(args='')

class TestCreate(TestBase):

    def test01(self):
        self.pkgname = self.outfile+'.refpkg'
        shutil.rmtree(self.pkgname, ignore_errors=True)
        self.cmd_ok(args='create --package-name=%(pkgname)s' % self)
        self.cmd_fails(args='create --package-name=%(pkgname)s' % self)

    def test02(self):
        self.taxonomy = os.path.join(datadir, 'tax_table.csv')
        self.seq_info = os.path.join(datadir, 'bv_refdata.csv')
        self.pkgname = self.outfile+'.refpkg'

        shutil.rmtree(self.pkgname, ignore_errors=True)
        self.cmd_ok(args='create --package-name=%(pkgname)s --taxonomy=%(taxonomy)s --seq-info=%(seq_info)s' % self)

    def test03(self):
        self.taxonomy = os.path.join(datadir, 'tax_table.csv')
        self.seq_info = os.path.join(datadir, 'bv_refdata.csv')
        self.tree_stats = os.path.join(datadir, 'phyml_aa_stats.txt')
        self.pkgname = self.outfile+'.refpkg'

        shutil.rmtree(self.pkgname, ignore_errors=True)
        self.cmd_ok(args='create --package-name=%(pkgname)s --taxonomy=%(taxonomy)s --seq-info=%(seq_info)s --tree-stats=%(tree_stats)s' % self)







