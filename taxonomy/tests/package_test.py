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

from Taxonomy.package import MLOutputParser

log = logging

module_name = os.path.split(sys.argv[0])[1].rstrip('.py')
outputdir = os.path.abspath(config.outputdir)
datadir = os.path.abspath(config.datadir)

class TestMLOutputParser(unittest.TestCase):

    def setUp(self):
        self.funcname = '_'.join(self.id().split('.')[-2:])

    # Verify we can extract data from all of our example input files.
    def test01(self):
        test_files = [os.path.join(datadir,'phyml_aa_stats.txt'),
                      os.path.join(datadir,'phyml_dna_stats.txt'),
                      os.path.join(datadir,'RAxML_info.re-estimated'),
                      os.path.join(datadir,'RAxML_info.testNuc')]
        for file_name in test_files:
            parser = MLOutputParser(file_name)
            result = parser.parse_ml_data()
            self.assertTrue(result)
