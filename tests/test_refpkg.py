import unittest
import tempfile
import shutil
import json
import copy
import os
import os.path

from taxtastic import refpkg
from . import config

class TestRefpkg(unittest.TestCase):
    maxDiff = None

    def test_fails_on_file_exists(self):
        # Trying to attach to a normal file as a Refpkg should fail.
        with tempfile.NamedTemporaryFile() as tf:
            self.assertRaises(ValueError, lambda: refpkg.Refpkg(tf.name))

    def test_fails_on_dir_without_manifest(self):
        # Trying to attach to a non-Refpkg directory should fail.
        with config.tempdir() as td:
            self.assertRaises(ValueError, lambda: refpkg.Refpkg(td))

    def test_create(self):
        # Attaching to an empty directory should create a new, empty Refpkg.
        scratch = tempfile.mkdtemp()
        try:
            pkg_path = os.path.join(scratch, 'test.refpkg')
            r = refpkg.Refpkg(pkg_path)
            self.assertEqual(r.contents, refpkg.manifest_template())
            self.assertEqual(os.listdir(pkg_path), ['CONTENTS.json'])
        finally:
            shutil.rmtree(scratch)

    def test_catch_invalid_json_keys(self):
        # Make sure a CONTENTS.json with invalid keys is caught
        scratch = tempfile.mkdtemp()
        try:
            pkg_path = os.path.join(scratch, 'test.refpkg')
            os.mkdir(pkg_path)
            with open(os.path.join(pkg_path, 'CONTENTS.json'), 'w') as h:
                json.dump({'boris': {}, 'hilda': {}}, h)
                self.assertRaises(ValueError,
                                  lambda: refpkg.Refpkg(pkg_path))
        finally:
            shutil.rmtree(scratch)

    def test_catch_nonmatching_md5(self):
        scratch = tempfile.mkdtemp()
        try:
            pkg_path = os.path.join(scratch, 'test.refpkg')
            os.mkdir(pkg_path)
            boris = os.path.join(pkg_path,'boris')
            with open(boris, 'w') as h:
                print >>h, "Hello, world!"
            with open(os.path.join(pkg_path, refpkg.Refpkg._manifest_name), 'w') as h:
                json.dump({'metadata': {},
                           'log': [],
                           'rollback': None,
                           'rollforward': None,
                           'files': {'meep': 'boris'},
                           'md5': {'meep': 0}},
                          h)
                self.assertRaises(ValueError, lambda: refpkg.Refpkg(pkg_path))
        finally:
            shutil.rmtree(scratch)

    def test_catch_nonmatching_md5_and_file_keys(self):
        scratch = tempfile.mkdtemp()
        try:
            pkg_path = os.path.join(scratch, 'test.refpkg')
            os.mkdir(pkg_path)
            boris = os.path.join(pkg_path,'boris')
            with open(boris, 'w') as h:
                print >>h, "Hello, world!"
            with open(os.path.join(pkg_path, refpkg.Refpkg._manifest_name), 'w') as h:
                json.dump({'metadata': {},
                           'log': [], 'rollback': None, 'rollforward': None,
                           'files': {'meep': 'boris', 'hilda': 'boris'},
                           'md5': {'meep': 0}},
                          h)
                self.assertRaises(ValueError, lambda: refpkg.Refpkg(pkg_path))
        finally:
            shutil.rmtree(scratch)

    def test_update_file(self):
        scratch = tempfile.mkdtemp()
        try:
            pkg_path = os.path.join(scratch, 'test.refpkg')
            r = refpkg.Refpkg(pkg_path)
            test_file = config.data_path('bv_refdata.csv')
            test_name = 'bv_refdata.csv'
            md5_value = refpkg.md5file(test_file)
            self.assertEqual(None, r.update_file('a', test_file))
            # Make sure it's properly written
            with open(os.path.join(r.path, r._manifest_name)) as h:
                self.assertEqual(json.load(h), r.contents)

            self.assertIn('a', r.contents['files'])
            self.assertEqual(r.file_name('a'), test_name)
            self.assertEqual(r.file_md5('a'), md5_value)

            self.assertEqual(None, r.update_file('b', test_file))
            self.assertEqual(r.file_name('b'), test_name + '1')
            self.assertEqual(r.file_md5('b'), md5_value)

            test_file2 = config.data_path('taxids1.txt')

            old_path = r.file_abspath('a')
            self.assertEqual(old_path, r.update_file('a', test_file2))
            self.assertTrue(os.path.exists(os.path.join(r.path, test_name)))
        finally:
            shutil.rmtree(scratch)

    def test_update_metadata(self):
        scratch = tempfile.mkdtemp()
        try:
            pkg_path = os.path.join(scratch, 'test.refpkg')
            r = refpkg.Refpkg(pkg_path)
            self.assertEqual(r.update_metadata('a', 'boris'), None)
            self.assertEqual(r.update_metadata('a', 'meep'), 'boris')
            with open(os.path.join(r.path, r._manifest_name)) as h:
                self.assertEqual(json.load(h), r.contents)
            self.assertEqual(r.metadata('a'), 'meep')
            self.assertEqual(r.metadata('b'), None)
        finally:
            shutil.rmtree(scratch)

    def test_reroot(self):
        with config.tempdir() as d:
            rpkg = os.path.join(d, 'reroot.refpkg')
            shutil.copytree(config.data_path('lactobacillus2-0.2.refpkg'), rpkg)
            r = refpkg.Refpkg(rpkg)
            r.reroot()
            self.assertEqual('9bdbf22f8bf140074d126f3d27989100',
                             r.file_md5('tree'))
            self.assertEqual(r.log(), ['Rerooting refpkg'])

    def test_pretend_reroot(self):
        with config.tempdir() as d:
            rpkg = os.path.join(d, 'reroot.refpkg')
            shutil.copytree(config.data_path('lactobacillus2-0.2.refpkg'), rpkg)
            r = refpkg.Refpkg(rpkg)
            r.reroot(pretend=True)
            self.assertEqual('2f11faa616fc7f04d7694436b5cca05f',
                             r.file_md5('tree'))

    def test_transaction(self):
        with config.tempdir() as d:
            rpkg = os.path.join(d, 'test.refpkg')
            shutil.copytree(config.data_path('lactobacillus2-0.2.refpkg'), rpkg)
            r = refpkg.Refpkg(rpkg)
            self.assertEqual(r.update_metadata('author', 'Boris and Hilda'),
                             "Noah Hoffman <ngh2@uw.edu>, Sujatha Srinivasan <ssriniva@fhcrc.org>, Erick Matsen <matsen@fhcrc.org>")
            self.assertEqual(r.current_transaction, None)
            self.assertEqual(r.log(),
                             ['Updated metadata: author=Boris and Hilda'])
            self.assertTrue(isinstance(r.contents['rollback'], dict))
            self.assertFalse('log' in r.contents['rollback'])

            original_log = copy.deepcopy(r.log())
            r.start_transaction()
            r.update_metadata('boris', 'meep')
            r.update_metadata('hilda', 'vrrp')
            r._log("Meep!")
            r.commit_transaction()
            self.assertFalse('boris' in r.contents['rollback']['metadata'])
            self.assertFalse('hilda' in r.contents['rollback']['metadata'])
            self.assertEqual(r.log(), ["Meep!"]+original_log)


    def test_failed_transaction(self):
        with config.tempdir() as d:
            rpkg = os.path.join(d, 'test.refpkg')
            shutil.copytree(config.data_path('lactobacillus2-0.2.refpkg'), rpkg)
            r = refpkg.Refpkg(rpkg)
            r = refpkg.Refpkg(rpkg)
            v = copy.deepcopy(r.contents)
            self.assertRaises(Exception,
                              lambda: r.update_file('tiddlywinks',
                                                    '/path/to/nonexistant/thing'))
            self.assertEqual(v, r.contents)

    def test_rollback(self):
        with config.tempdir() as d:
            rpkg = os.path.join(d, 'test.refpkg')
            shutil.copytree(config.data_path('lactobacillus2-0.2.refpkg'), rpkg)
            r = refpkg.Refpkg(rpkg)
            self.assertRaises(ValueError,
                              lambda: r.rollback())
            self.maxDiff = None
            v0 = copy.deepcopy(r.contents)
            self.assertFalse('boris' in r.contents['metadata'])
            r.start_transaction()
            r.update_metadata('boris', 'meep')
            r.update_file('boris', config.data_path('taxids1.txt'))
            r.commit_transaction()
            boris_path = r.file_abspath('boris')
            self.assertTrue('boris' in r.contents['files'])
            self.assertFalse('boris' in r.contents['rollback']['files'])

            self.assertFalse('boris' in r.contents['rollback']['metadata'])
            self.assertTrue('boris' in r.contents['metadata'])


            v1 = copy.deepcopy(r.contents)
            r.rollback()
            self.assertFalse('boris' in r.contents['files'])
            self.assertFalse('boris' in r.contents['md5'])
            self.assertTrue(os.path.exists(boris_path))
            v3 = copy.deepcopy(r.contents)
            v0.pop('rollforward')
            v3.pop('rollforward')
            self.assertEqual(v0, v3)
            r.rollforward()
            self.assertEqual(v1, r.contents)

            # We shouldn't be able to roll forward after running an unrelated operation
            r.rollback()
            r.update_metadata('boris', 'hilda')
            self.assertRaises(ValueError,
                              lambda: r.rollforward())

    def test_strip(self):
        with config.tempdir() as d:
            rpkg = os.path.join(d, 'test.refpkg')
            shutil.copytree(config.data_path('lactobacillus2-0.2.refpkg'), rpkg)
            r = refpkg.Refpkg(rpkg)
            self.assertFalse('boris' in r.contents['files'])
            r.update_file('boris', config.data_path('taxids1.txt'))
            boris_path = r.file_abspath('boris')
            r.rollback()
            self.assertFalse('boris' in r.contents['files'])
            original_log= r.log()
            r.strip()

            self.assertFalse('boris' in r.contents['files'])
            self.assertEqual(r.log(), ['Stripped refpkg (removed 1 files)'] + original_log)
            self.assertFalse(os.path.exists(boris_path))
            self.assertFalse(r.is_invalid())
            self.assertEqual(len(r.contents['files']), len(os.listdir(r.path))-1)

    def test_is_ill_formed(self):
        with config.tempdir() as d:
            rpkg = os.path.join(d, 'test.refpkg')
            shutil.copytree(config.data_path('lactobacillus2-0.2.refpkg'), rpkg)
            r = refpkg.Refpkg(rpkg)
            self.assertFalse(r.is_ill_formed())
            r.update_file('aln_fasta', config.data_path('little.fasta'))
            self.assertTrue(isinstance(r.is_ill_formed(), basestring))


if __name__ == '__main__':
    unittest.main()


