#!/usr/bin/env python

import sys
import os
import unittest
import logging
import pprint
import config
import collections
import json

from taxonomy.package import StatsParser

log = logging

module_name = os.path.split(sys.argv[0])[1].rstrip('.py')
outputdir = os.path.abspath(config.outputdir)
datadir = os.path.abspath(config.datadir)

test_files = ['phyml_aa_stats.txt',
              'phyml_dna_stats.txt',
              'RAxML_info.re-estimated',
              'RAxML_info.aa',
              'RAxML_info.testNuc']

test_paths = [os.path.join(datadir, f) for f in test_files]

class TestStatsParser(unittest.TestCase):

    def setUp(self):
        self.funcname = '_'.join(self.id().split('.')[-2:])

    def testRead(self):
        """
        Verify we can extract data from each input file.
        """

        for file_name in test_paths:
            log.info('reading %s' % file_name)
            parser = StatsParser(file_name)
            result = parser.parse_stats_data()
            self.assertTrue(result)

    def testGetValues(self):
        """
        Verify that each file results in values containing data.
        """

        for file_name in test_paths:
            log.info('reading %s' % file_name)
            parser = StatsParser(file_name)
            parser.parse_stats_data()
            values = parser.get_stats_values()
            self.assertTrue(isinstance(values, collections.defaultdict))
            # values shoud not be empty
            self.assertTrue(bool(dict(values)))

    def testWriteStatsJson(self):
        """
        Verify that JSON output can be written for each input file -
        no checks of content here.
        """
        for f in test_files:
            infile = os.path.join(datadir, f)
            outfile = os.path.join(outputdir, f)+'.json'
            try:
                os.remove(outfile)
            except OSError:
                pass

            parser = StatsParser(infile)
            parser.parse_stats_data()
            parser.write_stats_json(outfile)
            self.assertTrue(os.path.isfile(outfile))

    def testJsonRegression(self):
        """
        Regression testing of values parsed from each input file to
        its stored json equivalent.
        """

        def compare(d1, d2):
            """
            Recursively compare key, value dictionary pairs. Raises
            error on difference.
            """

            for k in set(d1.keys()) | set(d2.keys()):
                try:
                    v1, v2 = d1[k], d2[k]
                except KeyError, msg:
                    log.error(pprint.pformat(d1))
                    log.error(pprint.pformat(d2))
                    raise KeyError(msg)

                log.debug('%s %s %s' % (k, v1, v2))

                if isinstance(v1, dict):
                    compare(v1, v2)
                elif isinstance(v1, float):
                    self.assertAlmostEqual(v1, v2)
                else:
                    self.assertEqual(v1, v2)

        for statsfile in test_paths:
            log.info(statsfile)
            parser = StatsParser(statsfile)
            parser.parse_stats_data()
            svals = dict(parser.get_stats_values())
            with open(statsfile + '.json') as fobj:
                jvals = json.load(fobj)

            compare(svals, jvals)

