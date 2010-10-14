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

import config
import Taxonomy

log = logging

outputdir = os.path.abspath(config.outputdir)
datadir = os.path.abspath(config.datadir)

class TestBriansClass(unittest.TestCase):

    def setUp(self):
        self.funcname = '_'.join(self.id().split('.')[-2:])

    def tearDown(self):
        pass

    def test01(self):
        print datadir
        print outputdir
