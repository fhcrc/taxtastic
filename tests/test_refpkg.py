import unittest
import tempfile
import shutil
import json
import copy
import sys
import os

sys.path.insert(1, '../')
from taxtastic import refpkg

class TestRefpkg(unittest.TestCase):
    def test_fails_on_file_exists(self):
        # Trying to attach to a normal file as a Refpkg should fail.
        self.assertRaises(ValueError, lambda: refpkg.Refpkg('test_refpkg.py'))

    def test_fails_on_dir_without_manifest(self):
        # Trying to attach to a non-Refpkg directory should fail.
        self.assertRaises(ValueError, lambda: refpkg.Refpkg('../taxtastic'))

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
            test_file = '../testfiles/bv_refdata.csv'
            test_name = 'bv_refdata.csv'
            md5_value = refpkg.md5file(test_file)
            self.assertEqual(r.update_file('a', test_file),
                             ('a', md5_value))
            # Make sure it's properly written
            with open(os.path.join(r.path, r._manifest_name)) as h:
                self.assertEqual(json.load(h), r.contents)

            self.assertTrue('a' in r.contents['files'])
            self.assertEqual(r.file_name('a'), test_name)
            self.assertEqual(r.file_md5('a'), md5_value)

            self.assertEqual(r.update_file('b', test_file),
                             ('b', md5_value))
            self.assertEqual(r.file_name('b'), test_name + '1')
            self.assertEqual(r.file_md5('b'), md5_value)

            test_file2 = '../testfiles/taxids1.txt'
            md5_value2 = refpkg.md5file(test_file2)

            self.assertEqual(r.update_file('a', test_file2),
                             ('a', md5_value2))
            self.assertFalse(os.path.exists(os.path.join(r.path, test_name)))
            
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
        shutil.copytree('../testfiles/lactobacillus2-0.2.refpkg', '../testfiles/toreroot.refpkg')
        try:
            r = refpkg.Refpkg('../testfiles/toreroot.refpkg')
            r.reroot()
            self.assertEqual('cd6431ed636ddb12b83ea0e4f9713bee',
                             r.file_md5('tree'))
            self.assertEqual(r.log(), ['Rerooting refpkg'])
        finally:
            shutil.rmtree('../testfiles/toreroot.refpkg')

    def test_pretend_reroot(self):
        shutil.copytree('../testfiles/lactobacillus2-0.2.refpkg', '../testfiles/toreroot.refpkg')     
        try:
            r = refpkg.Refpkg('../testfiles/toreroot.refpkg')
            r.reroot(pretend=True)
            self.assertEqual('2f11faa616fc7f04d7694436b5cca05f',
                             r.file_md5('tree'))
        finally:
            shutil.rmtree('../testfiles/toreroot.refpkg')

    def test_transaction(self):
        shutil.copytree('../testfiles/lactobacillus2-0.2.refpkg', '../testfiles/test.refpkg')     
        try:
            r = refpkg.Refpkg('../testfiles/test.refpkg')
            self.assertEqual(r.update_metadata('author', 'Boris and Hilda'), 
                             "Noah Hoffman <ngh2@uw.edu>, Sujatha Srinivasan <ssriniva@fhcrc.org>, Erick Matsen <matsen@fhcrc.org>")
            self.assertEqual(r.current_transaction, None)
            self.assertEqual(r.log(), 
                             ['Updated metadata: author=Boris and Hilda'])
            self.assertTrue(isinstance(r.contents['rollback'], dict))
            self.assertFalse('log' in r.contents['rollback'])
        finally:
            shutil.rmtree('../testfiles/test.refpkg')

    def test_failed_transaction(self):
        shutil.copytree('../testfiles/lactobacillus2-0.2.refpkg', '../testfiles/test.refpkg')     
        try:
            r = refpkg.Refpkg('../testfiles/test.refpkg')
            v = copy.deepcopy(r.contents)
            self.assertRaises(Exception, 
                              lambda: r.update_file('tiddlywinks', 
                                                    '/path/to/nonexistant/thing'))
            self.assertEqual(v, r.contents)
        finally:
            shutil.rmtree('../testfiles/test.refpkg')

    def test_rollback(self):
        shutil.copytree('../testfiles/lactobacillus2-0.2.refpkg', 
                        '../testfiles/test.refpkg')
        try:
            r = refpkg.Refpkg('../testfiles/test.refpkg')
            self.assertRaises(ValueError,
                              lambda: r.rollback())
            self.maxDiff = None
            v0 = copy.deepcopy(r.contents)
            r.update_metadata('boris', 'meep')
            v1 = copy.deepcopy(r.contents)
            r.rollback()
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
        finally:
            shutil.rmtree('../testfiles/test.refpkg')
        


if __name__ == '__main__':
    unittest.main()

