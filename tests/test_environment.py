#!/usr/bin/env python

import os
import unittest
import logging

from . import config
import taxtastic

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
        version of the taxtastic package. Fail otherwise.
        """
        taxtastic_path = os.path.abspath(os.path.dirname(taxtastic.__file__))
        exp_path = os.path.abspath(os.path.join(outputdir, '..', 'taxtastic'))
        self.assertEqual(exp_path, taxtastic_path)
