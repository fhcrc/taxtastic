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

outputdir = path.abspath(config.outputdir)
datadir = path.abspath(config.datadir)

TestScriptBase.executable = './taxtest.py'
TestScriptBase.outputdir = outputdir
        
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
        
        # fails the second time because package already exists
        self.cmd_fails('create -P %(package)s -l 16s')
