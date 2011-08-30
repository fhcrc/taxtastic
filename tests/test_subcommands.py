import unittest
import tempfile
import shutil
import sys
import os

sys.path.insert(1, '../')
from taxtastic import refpkg
from taxtastic.subcommands import update, create

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
                metadata = False
            update.action(_Args())
            r._sync_from_disk()
            self.assertEqual(r.contents['files']['meep'], 'bv_refdata.csv')
            self.assertEqual(r.contents['files']['hilda'], 'bv_refdata.csv1')
        finally:
            shutil.rmtree(scratch)

    def test_metadata_action(self):
        scratch = tempfile.mkdtemp()
        try:
            pkg_path = os.path.join(scratch, 'test.refpkg')
            r = refpkg.Refpkg(pkg_path)
            class _Args(object):
                refpkg=pkg_path
                changes = ['meep=boris', 'hilda=vrrp']
                metadata = True
            update.action(_Args())
            r._sync_from_disk()
            self.assertEqual(r.metadata('meep'), 'boris')
            self.assertEqual(r.metadata('hilda'), 'vrrp')
        finally:
            shutil.rmtree(scratch)

class TestCreate(unittest.TestCase):
    def test_create(self):
        scratch = tempfile.mkdtemp()
        class _Args(object):
            clobber = True
            locus = 'Nowhere'
            description = 'A description'
            author = 'Boris the Mad Baboon'
            package_version = '0.3'
            package_name = os.path.join(scratch, 'test.refpkg')
            tree_stats = None
            aln_fasta = None
            aln_sto = None
            phylo_model = None
            seq_info = None
            mask = None
            profile = None
            readme = None
            tree = None
            taxonomy = None
        try:
            create.action(_Args())
            r = refpkg.Refpkg(_Args().package_name)
            self.assertEqual(r.metadata('locus'), 'Nowhere')
            self.assertEqual(r.metadata('description'), 'A description')
            self.assertEqual(r.metadata('author'), 'Boris the Mad Baboon')
            self.assertEqual(r.metadata('package_version'), '0.3')
            self.assertEqual(r.metadata('format_version'), '1.1')
        finally:
            shutil.rmtree(scratch)

if __name__ == '__main__':
    unittest.main()
