import contextlib
import unittest
import tempfile
import shutil
import copy
import sys
import os

from taxtastic import refpkg
from taxtastic.subcommands import update, create, strip, rollback, rollforward, taxtable, check

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
            self.assertEqual(r.contents['rollback'], None)

            args2 = _Args()
            args2.package_name = '/dev/null'
            args2.clobber = True
            self.assertEqual(1, create.action(args2))
        finally:
            shutil.rmtree(scratch)


class TestStrip(unittest.TestCase):
    def test_strip(self):
        shutil.copytree('../testfiles/lactobacillus2-0.2.refpkg', '../testfiles/tostrip.refpkg')
        try:
            r = refpkg.Refpkg('../testfiles/tostrip.refpkg')
            r.update_metadata('boris', 'hilda')
            r.update_metadata('meep', 'natasha')

            class _Args(object):
                refpkg = '../testfiles/tostrip.refpkg'
            strip.action(_Args())

            r._sync_from_disk()
            self.assertEqual(r.contents['rollback'], None)
            self.assertEqual(r.contents['rollforward'], None)
        finally:
            shutil.rmtree('../testfiles/tostrip.refpkg')

class TestRollback(unittest.TestCase):
    maxDiff = None
    def test_rollback(self):
        shutil.copytree('../testfiles/lactobacillus2-0.2.refpkg', '../testfiles/tostrip.refpkg')
        try:
            r = refpkg.Refpkg('../testfiles/tostrip.refpkg')
            original_contents = copy.deepcopy(r.contents)
            r.update_metadata('boris', 'hilda')
            r.update_metadata('meep', 'natasha')
            updated_contents = copy.deepcopy(r.contents)

            class _Args(object):
                refpkg = '../testfiles/tostrip.refpkg'
                def __init__(self, n):
                    self.n = n

            self.assertEqual(rollback.action(_Args(3)), 1)
            r._sync_from_disk()
            self.assertEqual(r.contents, updated_contents)

            self.assertEqual(rollback.action(_Args(2)), 0)
            r._sync_from_disk()
            self.assertEqual(r.contents['metadata'], original_contents['metadata'])
            self.assertEqual(r.contents['rollback'], None)
            self.assertNotEqual(r.contents['rollforward'], None)
        finally:
            shutil.rmtree('../testfiles/tostrip.refpkg')


class TestRollforward(unittest.TestCase):
    maxDiff = None
    def test_rollforward(self):
        shutil.copytree('../testfiles/lactobacillus2-0.2.refpkg', '../testfiles/tostrip.refpkg')
        try:
            r = refpkg.Refpkg('../testfiles/tostrip.refpkg')
            original_contents = copy.deepcopy(r.contents)
            r.update_metadata('boris', 'hilda')
            r.update_metadata('meep', 'natasha')
            updated_contents = copy.deepcopy(r.contents)
            r.rollback()
            r.rollback()

            class _Args(object):
                refpkg = '../testfiles/tostrip.refpkg'
                def __init__(self, n):
                    self.n = n

            self.assertEqual(rollforward.action(_Args(3)), 1)
            r._sync_from_disk()
            self.assertEqual(r.contents['metadata'], original_contents['metadata'])

            self.assertEqual(rollforward.action(_Args(2)), 0)
            r._sync_from_disk()
            self.assertEqual(r.contents['metadata'], updated_contents['metadata'])
            self.assertEqual(r.contents['rollforward'], None)
            self.assertNotEqual(r.contents['rollback'], None)
        finally:
            shutil.rmtree('../testfiles/tostrip.refpkg')


@contextlib.contextmanager
def scratch_file(unlink=True):
    """Create a temporary file and return its name.

    At the start of the with block a secure, temporary file is created
    and its name returned.  At the end of the with block it is
    deleted.
    """
    try:
        (tmp_fd, tmp_name) = tempfile.mkstemp(text=True)
        os.close(tmp_fd)
        yield tmp_name
    except ValueError, v:
        raise v
    else:
        if unlink:
            os.unlink(tmp_name)


class TestTaxtable(unittest.TestCase):
    maxDiff = None
    def test_taxids(self):
        with scratch_file() as out:
            with open(out, 'w') as h:
                class _Args(object):
                    database_file = '../testfiles/taxonomy.db'
                    taxids = '../testfiles/taxids1.txt'
                    taxnames = None
                    verbosity = 0
                    out_file = h
                self.assertEqual(taxtable.action(_Args()), 0)
            self.assertEqual(refpkg.md5file(out), '39d9b3c9250403bf2add6c5908f7c620')

    def test_invalid_taxid(self):
        with scratch_file() as out:
            with open(out, 'w') as h:
                class _Args(object):
                    database_file = '../testfiles/taxonomy.db'
                    taxids = 'horace,hilda'
                    taxnames = None
                    verbosity = 0
                    out_file = h
                self.assertEqual(taxtable.action(_Args()), 1)

class TestCheck(unittest.TestCase):
    def test_runs(self):
        class _Args(object):
            refpkg = '../testfiles/lactobacillus2-0.2.refpkg'
        self.assertEqual(check.action(_Args()), 0)

if __name__ == '__main__':
    unittest.main()
