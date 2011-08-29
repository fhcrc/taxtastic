import unittest
import tempfile
import shutil
import sys
import os

sys.path.insert(1, '../')
from taxtastic import refpkg
from taxtastic.subcommands import update

class TestUpdate(unittest.TestCase):
    def test_action(self):
        scratch = tempfile.mkdtemp()
        try:
            pkg_path = os.path.join(scratch, 'test.refpkg')
            r = refpkg.Refpkg(pkg_path)
            test_file = '../testfiles/bv_refdata.csv'
            class _Args(object):
                refpkg=pkg_path
                changes = ['meep='+test_file, 'hilda='+test_file]
            update.action(_Args())
            r.reread()
            self.assertEqual(r.contents['files']['meep'], 'bv_refdata.csv')
            self.assertEqual(r.contents['files']['hilda'], 'bv_refdata.csv1')
        finally:
            shutil.rmtree(scratch)

if __name__ == '__main__':
    unittest.main()
