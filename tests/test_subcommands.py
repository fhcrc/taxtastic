import sys; sys.path.insert(0, '../')
import contextlib
import unittest
import tempfile
import shutil
import copy
import os

from taxtastic import refpkg
from taxtastic.lonely import Tree
from taxtastic.subcommands import update, create, strip, rollback, rollforward, taxtable, check, lonelynodes, findcompany

import config
from config import OutputRedirectMixin

class TestUpdate(OutputRedirectMixin, unittest.TestCase):
    def test_action(self):
        with config.tempdir() as scratch:
            pkg_path = os.path.join(scratch, 'test.refpkg')
            r = refpkg.Refpkg(pkg_path)
            test_file = config.data_path('bv_refdata.csv')
            class _Args(object):
                refpkg=pkg_path
                changes = ['meep='+test_file, 'hilda='+test_file]
                metadata = False
            update.action(_Args())
            r._sync_from_disk()
            self.assertEqual(r.contents['files']['meep'], 'bv_refdata.csv')
            self.assertEqual(r.contents['files']['hilda'], 'bv_refdata.csv1')

    def test_metadata_action(self):
        with config.tempdir() as scratch:
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

class TestCreate(OutputRedirectMixin, unittest.TestCase):
    def test_create(self):
        with config.tempdir() as scratch:
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
            create.action(_Args())
            r = refpkg.Refpkg(_Args().package_name)
            self.assertEqual(r.metadata('locus'), 'Nowhere')
            self.assertEqual(r.metadata('description'), 'A description')
            self.assertEqual(r.metadata('author'), 'Boris the Mad Baboon')
            self.assertEqual(r.metadata('package_version'), '0.3')
            self.assertEqual(r.metadata('format_version'), '1.1')
            self.assertEqual(r.contents['rollback'], None)

            args2 = _Args()
            args2.package_name = os.path.join(scratch, 'test.refpkg')
            args2.clobber = True
            self.assertEqual(0, create.action(args2))


class TestStrip(OutputRedirectMixin, unittest.TestCase):
    def test_strip(self):
        with config.tempdir() as scratch:
            rpkg = os.path.join(scratch, 'tostrip.refpkg')
            shutil.copytree(config.data_path('lactobacillus2-0.2.refpkg'), rpkg)
            r = refpkg.Refpkg(rpkg)
            r.update_metadata('boris', 'hilda')
            r.update_metadata('meep', 'natasha')

            class _Args(object):
                refpkg = rpkg
            strip.action(_Args())

            r._sync_from_disk()
            self.assertEqual(r.contents['rollback'], None)
            self.assertEqual(r.contents['rollforward'], None)

class TestRollback(OutputRedirectMixin, unittest.TestCase):
    maxDiff = None
    def test_rollback(self):
        with config.tempdir() as scratch:
            rpkg = os.path.join(scratch, 'tostrip.refpkg')
            shutil.copytree(config.data_path('lactobacillus2-0.2.refpkg'), rpkg)
            r = refpkg.Refpkg(rpkg)
            original_contents = copy.deepcopy(r.contents)
            r.update_metadata('boris', 'hilda')
            r.update_metadata('meep', 'natasha')
            updated_contents = copy.deepcopy(r.contents)

            class _Args(object):
                refpkg = rpkg
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

class TestRollforward(OutputRedirectMixin, unittest.TestCase):
    maxDiff = None
    def test_rollforward(self):
        with config.tempdir() as scratch:
            rpkg = os.path.join(scratch, 'tostrip.refpkg')
            shutil.copytree(config.data_path('lactobacillus2-0.2.refpkg'), rpkg)
            r = refpkg.Refpkg(rpkg)
            original_contents = copy.deepcopy(r.contents)
            r.update_metadata('boris', 'hilda')
            r.update_metadata('meep', 'natasha')
            updated_contents = copy.deepcopy(r.contents)
            r.rollback()
            r.rollback()

            class _Args(object):
                refpkg = rpkg
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


class TestTaxtable(OutputRedirectMixin, unittest.TestCase):
    maxDiff = None
    @unittest.skip('Output varies.')
    def test_taxids(self):
        with scratch_file() as out:
            with open(out, 'w') as h:
                class _Args(object):
                    database_file = config.ncbi_master_db
                    taxids = config.data_path('taxids1.txt')
                    taxnames = None
                    seq_info = None
                    verbosity = 0
                    out_file = h
                self.assertEqual(taxtable.action(_Args()), 0)
            self.assertEqual(refpkg.md5file(out), '88ed5643d4754c60d5c472ff0b298f0f')

    def test_invalid_taxid(self):
        with scratch_file() as out:
            with open(out, 'w') as h:
                class _Args(object):
                    database_file = config.ncbi_master_db
                    taxids = 'horace,hilda'
                    taxnames = None
                    seq_info = None
                    verbosity = 0
                    out_file = h
                self.assertEqual(taxtable.action(_Args()), 1)

    def test_seqinfo(self):
        with tempfile.TemporaryFile() as tf, \
             open(config.data_path('simple_seqinfo.csv')) as ifp:
            class _Args(object):
                database_file = config.ncbi_master_db
                taxids = None
                taxnames = None
                seq_info = ifp
                out_file = tf
                verbosity = 0
            self.assertEqual(taxtable.action(_Args()), 0)
            # No output check at present
            self.assertTrue(tf.tell() > 0)

class TestCheck(OutputRedirectMixin, unittest.TestCase):
    def test_runs(self):
        class _Args(object):
            refpkg = config.data_path('lactobacillus2-0.2.refpkg')
        self.assertEqual(check.action(_Args()), 0)

def test_lonelynodes(capsys,tmpdir):
    t = Tree(1, rank='phylum', tax_name='a')(
             Tree(3, rank='order', tax_name='b')(
                 Tree(4, rank='class', tax_name='c')(
                     Tree(5, rank='family', tax_name='d')(
                         Tree(7, rank='genus', tax_name='e')),
                     Tree(6, rank='family', tax_name='f'))))
    infile = tmpdir.join('junk.taxtable')
    expected = ''.join(["3 # order b\n",
                          "4 # class c\n",
                          "7 # genus e\n"])
    taxtable = '''"tax_id","parent_id","rank","tax_name"
"1","1","phylum","a"
"3","1","order","b"
"4","3","class","c"
"5","4","family","d"
"7","5","genus","e"
"6","4","family","f"'''
    with open(str(infile), 'w') as h:
            print >>h, taxtable
    class _Args(object):
        target = str(infile)
        output = None
        verbose = True
    status = lonelynodes.action(_Args())
    assert status == 0
    out, err = capsys.readouterr()
    assert out == expected

def test_findcompany(capsys):
    class _Args(object):
        taxdb = '../testfiles/small_taxonomy.db'
        tax_ids = ['1239', '186801']
        input = None
        output = None
        cut = True
    status = findcompany.action(_Args())
    out, err = capsys.readouterr()
    assert status == 0
    assert err == ""
    assert out.strip() == "562\n1280"




