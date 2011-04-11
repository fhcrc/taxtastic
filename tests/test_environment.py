#!/usr/bin/env python

import os
import unittest
import logging

import config
import taxonomy

log = logging

outputdir = os.path.abspath(config.outputdir)
datadir = os.path.abspath(config.datadir)

class TestWhereWeAre(unittest.TestCase):

    def setUp(self):
        self.funcname = '_'.join(self.id().split('.')[-2:])

    def tearDown(self):
        pass

    def testPackageLocation(self):
        """
        We're assuming that unit tests are being run using local
        version of the taxonomy package. Fail otherwise.
        """
        
        self.assertTrue(taxonomy.__file__.startswith('taxonomy/__init__.py'))
                
