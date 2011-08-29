import itertools
import hashlib
import shutil
import os
import sqlite3
import json
import csv

def md5file(path):
    md5 = hashlib.md5()
    with open(path) as h:
        for block in iter(lambda: h.read(4096), ''):
            md5.update(block)
    return md5.hexdigest()


class VerificationFailed(Exception):
    pass

class NoAncestor(Exception):
    pass

class OnUpdate(object):
    def __init__(self, proxied):
        self.proxied = proxied
        self.setter = None

    def __get__(self, inst, cls):
        if inst is None:
            return self
        return getattr(inst, self.proxied)

    def __set__(self, inst, value):
        if self.setter:
            self.setter(inst, value)
        setattr(inst, self.proxied, value)

    def __call__(self, f):
        self.setter = f
        return self

class _IntermediateTaxon(object):
    def __init__(self, tax_id, parent, rank, tax_name):
        self.children = set()
        self.tax_id = tax_id
        self.parent = parent
        self.rank = rank
        self.tax_name = tax_name

    _parent = _adjacent_to = None

    @OnUpdate('_parent')
    def parent(self, p):
        if self.parent is not None:
            self.parent.children.discard(self)
        if p is not None:
            p.children.add(self)

    @OnUpdate('_adjacent_to')
    def adjacent_to(self, n):
        if n is not None:
            self.parent = n.parent

    def iterate_children(self, on_pop=None, including_self=True):
        search_stack = [(None, set([self]))]
        while search_stack:
            if not search_stack[-1][1]:
                parent, _ = search_stack.pop()
                if on_pop:
                    on_pop(parent)
                continue
            node = search_stack[-1][1].pop()
            if node is not self or including_self:
                yield node
            search_stack.append((node, node.children.copy()))

class Refpkg(object):
    _manifest_name = 'CONTENTS.json'
    _manifest_template = {'metadata': {},
                        'files': {},
                        'md5': {}}

    def __init__(self, path):
        """Create a reference to a new or existing RefPkg at *path*.

        If there is already a RefPkg at *path*, a reference is
        returned to that RefPkg.  If *path* does not exist, then an
        empty RefPkg is created.
        """
        # The logic of __init__ is complicated by having to check for
        # validity of a refpkg.  Much of its can be dispatched to the
        # isvalid method, but I want that to work at any time on the
        # RefPkg object, so it must have the RefPkg's manifest already
        # in place.
        self.path = os.path.abspath(path)
        if not(os.path.exists(path)):
            os.mkdir(path)
            with open(os.path.join(path, self._manifest_name), 'w') as h:
                json.dump(self._manifest_template, h)
        if not(os.path.isdir(path)):
            raise ValueError("%s is not a valid RefPkg" % (path,))
        # Try to load the Refpkg and check that it's valid
        manifest_path = os.path.join(path, self._manifest_name)
        if not(os.path.isfile(manifest_path)):
            raise ValueError(("%s is not a valid RefPkg - "
                              "could not find manifest file %s") % \
                                 (path, self._manifest_name))
        with open(manifest_path) as h:
            self.contents = json.load(h)
        error = self.isinvalid()
        if error:
            raise ValueError("%s is not a valid RefPkg: %s" % (path, error))


    def isinvalid(self):
        """Check if this RefPkg is invalid.

        Valid means that it contains a properly named manifest, and
        each of the files described in the manifest exists and has the
        proper MD5 hashsum.

        If the Refpkg is valid, isinvalid returns False.  Otherwise it
        returns a nonempty string describing the error.
        """
        # Contains a manifest file
        if not(os.path.isfile(os.path.join(self.path, self._manifest_name))):
            return "No manifest file %s found" % self._manifest_name
        # Manifest file contains the proper keys
        for k in self._manifest_template.keys():
            if not(k in self.contents):
                return "Manifest file missing key %s" % k
            if not(isinstance(self.contents[k], dict)):
                return "Key %s in manifest did not refer to a dictionary" % k
        # MD5 keys and filenames are in one to one correspondence
        if self.contents['files'].keys() != self.contents['md5'].keys():
            return ("Files and MD5 sums in manifest do not "
                    "match (files: %s, MD5 sums: %s)") % \
                    (self.contents['files'].keys(), 
                     self.contents['md5'].keys())
        # All files in the manifest exist and match the MD5 sums
        for key,filename in self.contents['files'].iteritems():
            expected_md5 = self.contents['md5'][key]
            filepath = os.path.join(self.path, filename)
            if not(os.path.exists(filepath)):
                return "File %s referred to by key %s not found in refpkg" % \
                    (filename, key)
            found_md5 = md5file(filepath)
            if found_md5 != expected_md5:
                return ("File %s referred to by key %s did "
                        "not match its MD5 sum (found: %s, expected %s)") % \
                        (found_md5, expected_md5)
        return False

    def _sync_to_disk(self):
        """Write any changes made on Refpkg to disk.

        Other methods of Refpkg that alter the contents of the package
        will call this method themselves.  Generally you should never
        have to call it by hand.  The only exception would be if
        another program has changed the Refpkg on disk while your
        program is running and you want to force your version over it.
        Otherwise it should only be called by other methods of refpkg.
        """
        with open(os.path.join(self.path, self._manifest_name), 'w') as h:
            json.dump(self.contents, h)

    def _sync_from_disk(self):
        """Read any changes made on disk to this Refpkg.

        This is necessary if other programs are making changes to the
        Refpkg on disk and your program must be synchronized to them.
        """
        with open(os.path.join(self.path, self._manifest_name)) as h:
            self.contents = json.load(h)
        error = self.isinvalid()
        if error:
            raise ValueError("Refpkg is invalid: %s" % error)

    def metadata(self, key):
        return self.contents['metadata'].get(key)

    def update_metadata(self, key, value):
        """Set *key* in the metadata to *value*.

        Returns the previous value of *key*, or None if the key was
        not previously set.
        """
        old_value = self.contents['metadata'].get(key)
        self.contents['metadata'][key] = value
        self._sync_to_disk()
        return old_value

    def update_file(self, key, new_path):
        """Insert file *new_path* into the Refpkg under *key*.

        The filename of *new_path* will be preserved in the Refpkg
        unless it would conflict with a previously existing file, in
        which case a suffix is appended which makes it unique.  Any
        file previously referred to by *key* is deleted.
        """
        if not(os.path.isfile(new_path)):
            raise ValueError("Cannot update Refpkg with file %s" % (new_path,))
        md5_value = md5file(new_path)
        filename = os.path.basename(new_path)
        while os.path.exists(os.path.join(self.path, filename)):
            filename += "1"
        if key in self.contents['files']:
            os.unlink(os.path.join(self.path, self.contents['files'][key]))
        shutil.copyfile(new_path, os.path.join(self.path, filename))
        self.contents['files'][key] = filename
        self.contents['md5'][key] = md5_value
        self._sync_to_disk()
        return (key, md5_value)

    def file_abspath(self, key):
        """Return the absolute path to the file referenced by *key*."""
        return os.path.join(self.path, self.file_name(key))

    def file_name(self, key):
        """Return the name of the file referenced by *key* in the refpkg."""
        if not(key in self.contents['files']):
            raise ValueError("No such resource key %s in refpkg" % key)
        return self.contents['files'][key]

    def file_md5(self, key):
        """Return the MD5 sum of the file reference by *key*."""
        if not(key in self.contents['md5']):
            raise ValueError("No such resource key %s in refpkg" % key)
        return self.contents['md5'][key]


