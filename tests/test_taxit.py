#!/usr/bin/env python

import sys
import os
from os import path
import unittest
import logging
import shutil
import commands
import re

import config
from config import TestScriptBase, rmdir

log = logging

datadir = config.datadir

TestScriptBase.executable = './taxit'
TestScriptBase.outputdir = config.outputdir
TestScriptBase.ncbi_master_db = config.ncbi_master_db

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

    def test01(self):
        self.cmd_fails('create -P %(package)s')

    def test02(self):
        """
        Create a minimal package.
        """
        
        rmdir(self.package)
        self.cmd_ok('create -P %(package)s -l 16s')
        self.assertTrue(path.exists(self.package))
        self.assertTrue(path.exists(path.join(self.package,'CONTENTS.json')))

        # test the --clobber option
        self.cmd_ok('create -P %(package)s -l 16s --clobber')
        
        # fails without --clobber because package already exists
        self.cmd_fails('create -P %(package)s -l 16s')


class TestTaxTable(TestScriptBase):
    """
    Unit tests for the taxtable sub-command.
    """
    # def test01(self):
    #     """
    #     Minimal test: a existing database is opened and that's it.
    #     """
    #     self.cmd_ok('taxtable --database-file %(ncbi_master_db)s')

    def test02(self):
        """
        Invalid arguments should cause a failure.
        """
        self.cmd_fails('taxtable --not-an-argument')
