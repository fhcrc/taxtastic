"""Test cases for RAxml parsing code in taxtastic/utils.py."""
import os
import sys
import json
import unittest

#import config
sys.path.insert(1, '../')
from taxtastic.utils import *

datadir = os.path.abspath('../testfiles')

class TestRaxmlStatsParse(unittest.TestCase):
    def test_try_set_fields(self):
        self.assertEqual({}, try_set_fields({}, r'asdf', ''))
        self.assertEqual({}, try_set_fields({}, r'asdf', 'asdf'))
        self.assertEqual({'d': 'asdf'}, 
                         try_set_fields({}, r'(?P<d>asdf)', 'boris asdf him'))
        self.assertEqual({'a': 'asdf', 'b': 'hilda'},
                         try_set_fields({}, r'(?P<a>asdf) then (?P<b>hilda)', 
                                        'asdf then hilda'))

    def test_aa_parse(self):
        target = {
            "empirical_frequencies": True,
            "datatype": "AA", 
            "subs_model": "WAG", 
            "program": "RAxML version 7.2.6", 
            "ras_model": "gamma", 
            "gamma": {
                "alpha": 0.73156900000000002, 
                "n_cats": 4
                }
            }
        with open(os.path.join(datadir, 'RAxML_info.aa')) as f:
            found = parse_raxml(f)
        self.assertEqual(target, found)

    def test_re_estimated_raxml(self):
        target = {
            "empirical_frequencies": True,
            "datatype": "DNA", 
            "subs_model": "GTR", 
            "program": "RAxML version 7.2.6", 
            "ras_model": "gamma", 
            "gamma": {
                "alpha": 0.54468399999999995, 
                "n_cats": 4
                }, 
            "subs_rates": {
                "ac": 0.723082, 
                "gt": 1.0, 
                "ag": 1.779172, 
                "cg": 0.81035199999999996, 
                "at": 1.2083870000000001, 
                "ct": 3.4005359999999998
                }
            }
        with open(os.path.join(datadir, 'RAxML_info.re-estimated')) as f:
            found = parse_raxml(f)
        self.assertEqual(target, found)

    def test_nuc(self):
        target = {
            "empirical_frequencies": True,
            "datatype": "DNA", 
            "subs_model": "GTR", 
            "program": "RAxML version 7.2.7", 
            "ras_model": "gamma", 
            "gamma": {
                "alpha": 0.38467699999999999, 
                "n_cats": 4
                }, 
            "subs_rates": {
                "ac": 0.49105900000000002, 
                "gt": 1.0, 
                "ag": 1.2250099999999999, 
                "cg": 0.53913699999999998, 
                "at": 0.91300800000000004, 
                "ct": 2.6734360000000001
                }
            }
        with open(os.path.join(datadir, 'RAxML_info.testNuc')) as f:
            found = parse_raxml(f)
        self.assertEqual(target, found)

