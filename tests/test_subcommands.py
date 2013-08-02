import contextlib
import unittest
import tempfile
import shutil
import copy
import os
import os.path

from taxtastic import refpkg
from taxtastic.subcommands import update, create, strip, rollback, rollforward, taxtable, check, add_to_taxtable

import config
from config import OutputRedirectMixin

class TestUpdate(OutputRedirectMixin, unittest.TestCase):
    def test_action(self):
        with config.tempdir() as scratch:
            pkg_path = os.path.join(scratch, 'test.refpkg')
            r = refpkg.Refpkg(pkg_path, create=True)
            test_file = config.data_path('bv_refdata.csv')
            class _Args(object):
                refpkg=pkg_path
                changes = ['meep='+test_file, 'hilda='+test_file]
                metadata = False
            update.action(_Args())
            r._sync_from_disk()
            self.assertEqual(r.contents['files']['meep'], 'bv_refdata.csv')

            # Second file should have been assigned a non-clashing name
            h = r.contents['files']['hilda']
            self.assertNotEqual(h, 'bv_refdata.csv')
            self.assertTrue(h.startswith('bv_refdata'))
            self.assertTrue(h.endswith('.csv'))

            self.assertTrue(os.path.exists(r.resource_path('hilda')))

    def test_metadata_action(self):
        with config.tempdir() as scratch:
            pkg_path = os.path.join(scratch, 'test.refpkg')
            r = refpkg.Refpkg(pkg_path, create=True)
            class _Args(object):
                refpkg=pkg_path
                changes = ['meep=boris', 'hilda=vrrp']
                metadata = True
            update.action(_Args())
            r._sync_from_disk()
            self.assertEqual(r.metadata('meep'), 'boris')
            self.assertEqual(r.metadata('hilda'), 'vrrp')

class TestCreate(OutputRedirectMixin, unittest.TestCase):

    def setUp(self):
        super(TestCreate, self).setUp()
        class _Args(object):
            clobber = True
            locus = 'Nowhere'
            description = 'A description'
            author = 'Boris the Mad Baboon'
            package_version = '0.3'
            tree_stats = None
            stats_type = None
            aln_fasta = None
            aln_sto = None
            phylo_model = None
            seq_info = None
            mask = None
            profile = None
            readme = None
            tree = None
            taxonomy = None
            reroot = False
            rppr = 'rppr'
            def __init__(self, scratch):
                self.package_name = os.path.join(scratch, 'test.refpkg')
        self._Args = _Args

    def test_create(self):
        with config.tempdir() as scratch:
            args = self._Args(scratch)
            create.action(args)
            r = refpkg.Refpkg(args.package_name, create=False)
            self.assertEqual(r.metadata('locus'), 'Nowhere')
            self.assertEqual(r.metadata('description'), 'A description')
            self.assertEqual(r.metadata('author'), 'Boris the Mad Baboon')
            self.assertEqual(r.metadata('package_version'), '0.3')
            self.assertEqual(r.metadata('format_version'), '1.1')
            self.assertEqual(r.contents['rollback'], None)
            args2 = self._Args(scratch)
            args2.package_name = os.path.join(scratch, 'test.refpkg')
            args2.clobber = True
            self.assertEqual(0, create.action(args2))

    def _test_create_phylo_model(self, stats_path, stats_type=None):
        with config.tempdir() as scratch:
            args = self._Args(scratch)
            args.tree_stats = stats_path
            args.stats_type = stats_type
            create.action(args)

            r = refpkg.Refpkg(args.package_name, create=False)
            self.assertIn('phylo_model', r.contents['files'])

    def test_create_phyml_aa(self):
        stats_path = os.path.join(config.datadir, 'phyml_aa_stats.txt')
        self._test_create_phylo_model(stats_path)
        self._test_create_phylo_model(stats_path, 'PhyML')
        self.assertRaises(ValueError, self._test_create_phylo_model,
                          stats_path, 'FastTree')
        self.assertRaises(ValueError, self._test_create_phylo_model,
                          stats_path, 'garli')


class TestStrip(OutputRedirectMixin, unittest.TestCase):
    def test_strip(self):
        with config.tempdir() as scratch:
            rpkg = os.path.join(scratch, 'tostrip.refpkg')
            shutil.copytree(config.data_path('lactobacillus2-0.2.refpkg'), rpkg)
            r = refpkg.Refpkg(rpkg, create=False)
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
            r = refpkg.Refpkg(rpkg, create=False)
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
            r = refpkg.Refpkg(rpkg, create=False)
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

class TestAddToTaxtable(OutputRedirectMixin, unittest.TestCase):
    maxDiff = None

    def test_seqinfo(self):
        with tempfile.TemporaryFile() as tf, \
             open(config.data_path('minimal_taxonomy.csv')) as taxtable_fp, \
             open(config.data_path('minimal_add_taxonomy.csv')) as extra_nodes_fp:
            class _Args(object):
                extra_nodes_csv = extra_nodes_fp
                taxtable = taxtable_fp
                out_file = tf
                verbosity = 0
            self.assertFalse(add_to_taxtable.action(_Args()))
            # No output check at present
            self.assertTrue(tf.tell() > 0)

class TestCheck(OutputRedirectMixin, unittest.TestCase):
    def test_runs(self):
        class _Args(object):
            refpkg = config.data_path('lactobacillus2-0.2.refpkg')
        self.assertEqual(check.action(_Args()), 0)
