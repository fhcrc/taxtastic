import unittest
import tempfile
import shutil
import json
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
            self.assertEqual(r.contents, refpkg.Refpkg._manifest_template)
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


if __name__ == '__main__':
    unittest.main()

