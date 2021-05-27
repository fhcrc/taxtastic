"""
taxtastic/refpkg.py

Implements an object, Refpkg, for the creation and manipulation of
reference packages for pplacer.

Note that Refpkg objects are *NOT* thread safe!
"""
# This file is part of taxtastic.
#
#    taxtastic is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    taxtastic is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with taxtastic.  If not, see <http://www.gnu.org/licenses/>.
import errno
import os
import time
import warnings
import sys
import contextlib
import copy
import csv
import functools
import hashlib
import json
import shutil
import subprocess
import tempfile
import zipfile

import dendropy
from decorator import decorator
from fastalite import fastalite

from taxtastic import utils, taxdb

FORMAT_VERSION = '1.1'


class DerivedFileNotUpdatedWarning(UserWarning):
    pass


def is_string(val):
    if sys.version_info.major == 2:
        return isinstance(val, basestring)
    else:
        return isinstance(val, str)


def md5file(fobj):
    md5 = hashlib.md5()
    md5.update(fobj.read())
    return md5.hexdigest()


@contextlib.contextmanager
def scratch_file(unlink=True, **kwargs):
    """Create a temporary file and return its name.

    Additional arguments are passed to :class:`tempfile.NamedTemporaryFile`

    At the start of the with block a secure, temporary file is created
    and its name returned.  At the end of the with block it is
    deleted.
    """
    kwargs['delete'] = False
    tf = tempfile.NamedTemporaryFile(**kwargs)
    tf.close()
    try:
        yield tf.name
    finally:
        if unlink:
            os.unlink(tf.name)


def manifest_template():
    return {'metadata': {'create_date': time.strftime('%Y-%m-%d %H:%M:%S'),
                         'format_version': FORMAT_VERSION},
            'files': {},
            'md5': {},
            'log': [],
            'rollback': None,
            'rollforward': None}


# The transaction and rollback/rollforward system used by Refpkg uses
# a data structure from purely functional programming called a zipper.
# The current state is augmented with an ordered list of previous
# states and an orderd list of subsequent states.  Rolling back to a
# previous state of the refpkg means pushing the current state onto
# the ordered list of subsequent states, and popping the first of the
# previous states to become the new current state.  Rolling forward
# again runs in just the opposite direction.

# The log is maintained only on the current state to save space.  This
# slightly complicates Refpkg.rollback and Refpkg.rollforward.  Log
# messages for rollforward transactions are stored with the future
# states to make keeping the log up to date simple.

# In order to have a sensible set of states, we want compound commands
# to only produce a single state.  That is handled by this transaction
# decorator.  Since Refpkgs are already not threadsafe, we can use the
# Refpkg itself to hold the current transaction.  If there is no
# transaction, then we start one, and only that outermost transaction
# will write a log message.  Further, that outer most command will be
# reverted as a whole.  Logging is set by calling self._log(msg), but
# it must be the last command run which might contain transaction
# decorated calls!  Thus
#
#  self.update_metadata(...)
#  self._log(...)
#
# works properly, but
#  self._log(...)
#  self.update_metadata(...)
#
# does not.  The update_metadata call will overwrite the log that was
# set by _log.

@decorator
def transaction(f, self, *args, **kwargs):
    if self.current_transaction:
        return f(self, *args, **kwargs)
    else:
        self.start_transaction()
        try:
            r = f(self, *args, **kwargs)
            self.commit_transaction()
            return r
        except Exception:
            self.contents = self.current_transaction['rollback']
            self._sync_to_disk()
            raise
        finally:
            self.current_transaction = None


class NoAncestor(Exception):
    pass


class Refpkg(object):
    _manifest_name = 'CONTENTS.json'

    def __init__(self, path, create=None):
        """Create a reference to a new or existing RefPkg at *path*.

        If there is already a RefPkg at *path*, a reference is returned to that
        RefPkg. If *path* does not exist and *create* is true, then an empty
        RefPkg is created.
        """
        # The logic of __init__ is complicated by having to check for
        # validity of a refpkg.  Much of its can be dispatched to the
        # isvalid method, but I want that to work at any time on the
        # RefPkg object, so it must have the RefPkg's manifest already
        # in place.
        if create is None:
            warnings.warn(
                "the default value of True for the `create` parameter is "
                "deprecated; it will eventually become False",
                DeprecationWarning, stacklevel=2)
            create = True

        self.current_transaction = None
        self.path = os.path.abspath(path)
        if not os.path.exists(path):
            if create:
                os.mkdir(path)
                with self.open_manifest('w') as h:
                    json.dump(manifest_template(), h, indent=4)
                    h.write('\n')
            else:
                raise ValueError(
                    "Reference package {0} does not exist.".format(path))

        if zipfile.is_zipfile(path):
            self._install_zipfile_handlers()

        self._sync_from_disk()
        self._set_defaults()

        self.db = None

    def _install_zipfile_handlers(self):
        self._archive = zipfile.ZipFile(self.path)
        archive_dirs = [zi for zi in self._archive.infolist()
                        if zi.filename.endswith('/')]
        if len(archive_dirs) != 1:
            raise ValueError(
                'zipped reference packages must contain exactly one directory')
        archive_dir = archive_dirs[0].filename

        def zipopen(name, *mode):
            return self._archive.open(archive_dir + name, *mode)
        self.open = zipopen

        def file_path(name):
            with scratch_file(unlink=False) as path:
                src = self.open(name)
                with open(path, 'wb') as dst:
                    shutil.copyfileobj(src, dst)
            return path
        self.file_path = file_path

    # can be stubbed out to provide an alternative mechanism for
    # providing data (e.g. for testing)
    file_factory = open

    def file_path(self, name):
        """
        Return the full path to a file in the reference package.

        Do not use this method if it is at all possible to use the *open*
        method instead.
        """
        return os.path.join(self.path, name)

    def open(self, name, *mode):
        """
        Return an open file object for a file in the reference package.
        """
        return self.file_factory(self.file_path(name), *mode)

    def open_manifest(self, *mode):
        """
        Return an open file object for the manifest file in this reference
        package.
        """
        return self.open(self._manifest_name, *mode)

    def open_resource(self, resource, *mode):
        """
        Return an open file object for a particular named resource in this
        reference package.
        """
        return self.open(self.resource_name(resource), *mode)

    def resource_name(self, resource):
        """
        Return the name of the file within the reference package for a
        particular named resource.
        """
        if not(resource in self.contents['files']):
            raise ValueError("No such resource %r in refpkg" % (resource,))
        return self.contents['files'][resource]

    def resource_md5(self, resource):
        """Return the stored MD5 sum for a particular named resource."""
        if not(resource in self.contents['md5']):
            raise ValueError("No such resource %r in refpkg" % (resource,))
        return self.contents['md5'][resource]

    def calculate_resource_md5(self, resource):
        """Calculate the MD5 sum for a particular named resource."""
        with self.open_resource(resource, 'rb') as f:
            return md5file(f)

    def resource_path(self, resource):
        """
        Return the path to the file within the reference package for a
        particular named resource.

        Do not use this method if it is at all possible to use the
        *open_resource* method instead.
        """
        return self.file_path(self.resource_name(resource))

    def _set_defaults(self):
        """
        Set some default values in the manifest.

        This method should be called after loading from disk, but before
        checking the integrity of the reference package.
        """
        self.contents.setdefault('log', [])
        self.contents.setdefault('rollback', None)
        self.contents.setdefault('rollforward', None)

    def _sync_to_disk(self):
        """Write any changes made on Refpkg to disk.

        Other methods of Refpkg that alter the contents of the package
        will call this method themselves.  Generally you should never
        have to call it by hand.  The only exception would be if
        another program has changed the Refpkg on disk while your
        program is running and you want to force your version over it.
        Otherwise it should only be called by other methods of refpkg.
        """
        with self.open_manifest('w') as h:
            json.dump(self.contents, h, indent=4)
            h.write('\n')

    def _sync_from_disk(self):
        """Read any changes made on disk to this Refpkg.

        This is necessary if other programs are making changes to the
        Refpkg on disk and your program must be synchronized to them.
        """

        try:
            fobj = self.open_manifest('r')
        except IOError as e:
            if e.errno == errno.ENOENT:
                raise ValueError(
                    "couldn't find manifest file in %s" % (self.path,))
            elif e.errno == errno.ENOTDIR:
                raise ValueError("%s is not a directory" % (self.path,))
            else:
                raise
        with fobj:
            self.contents = json.load(fobj)

        self._set_defaults()
        self._check_refpkg()

    def _add_file(self, key, path):
        """Copy a file into the reference package."""
        filename = os.path.basename(path)
        base, ext = os.path.splitext(filename)
        if os.path.exists(self.file_path(filename)):
            with tempfile.NamedTemporaryFile(
                    dir=self.path, prefix=base, suffix=ext) as tf:
                filename = os.path.basename(tf.name)
        shutil.copyfile(path, self.file_path(filename))
        self.contents['files'][key] = filename

    def _delete_file(self, path):
        """Delete a file from the reference package."""
        os.unlink(self.file_path(path))

    def metadata(self, key):
        """Return the metadata value associated to *key*."""
        return self.contents['metadata'].get(key)

    def file_keys(self):
        """Return a list of all the keys referring to files in this refpkg."""
        return list(self.contents['files'].keys())

    def metadata_keys(self):
        """Return a list of all the keys referring to metadata in this refpkg."""
        return list(self.contents['metadata'].keys())

    def _log(self, msg):
        """Set the log message for this operation.

        When writing transactions that encapsulate several others,
        pass the log message to the call to ``commit_transaction``
        instead.
        """
        self.current_transaction['log'] = msg

    def log(self):
        """Returns the log of this refpkg.

        The log is a list of strings, one per operation, from newest
        to oldest.
        """
        return self.contents['log']

    def is_invalid(self):
        """Check if this RefPkg is invalid.

        Valid means that it contains a properly named manifest, and
        each of the files described in the manifest exists and has the
        proper MD5 hashsum.

        If the Refpkg is valid, is_invalid returns False.  Otherwise it
        returns a nonempty string describing the error.
        """
        # Manifest file contains the proper keys
        for k in ['metadata', 'files', 'md5']:
            if not(k in self.contents):
                return "Manifest file missing key %s" % k
            if not(isinstance(self.contents[k], dict)):
                return "Key %s in manifest did not refer to a dictionary" % k

        if not('rollback' in self.contents):
            return "Manifest file missing key rollback"
        if not(isinstance(self.contents['rollback'], dict)) and self.contents[
                "rollback"] is not None:
            return ("Key rollback in manifest did not refer to a "
                    "dictionary or None, found %s") % str(self.contents['rollback'])

        if not('rollforward' in self.contents):
            return "Manifest file missing key rollforward"

        if self.contents['rollforward'] is not None:
            if not(isinstance(self.contents['rollforward'], list)):
                return "Key rollforward was not a list, found %s" % str(
                    self.contents['rollforward'])
            elif len(self.contents['rollforward']) != 2:
                return "Key rollforward had wrong length, found %d" % \
                    len(self.contents['rollforward'])
            elif not is_string(self.contents['rollforward'][0]):
                print(type(self.contents['rollforward'][0]))
                return "Key rollforward's first entry was not a string, found %s" % \
                    str(self.contents['rollforward'][0])
            elif not(isinstance(self.contents['rollforward'][1], dict)):
                return "Key rollforward's second entry was not a dict, found %s" % \
                    str(self.contents['rollforward'][1])

        if not("log" in self.contents):
            return "Manifest file missing key 'log'"

        if not(isinstance(self.contents['log'], list)):
            return "Key 'log' in manifest did not refer to a list"

        # MD5 keys and filenames are in one to one correspondence
        if self.contents['files'].keys() != self.contents[
                'md5'].keys():
            return ("Files and MD5 sums in manifest do not "
                    "match (files: %s, MD5 sums: %s)") % \
                (list(self.contents['files'].keys()),
                 list(self.contents['md5'].keys()))

        # All files in the manifest exist and match the MD5 sums
        for key, filename in self.contents['files'].items():
            # we don't need to explicitly check for existence;
            # calculate_resource_md5 will open the file for us.
            expected_md5 = self.resource_md5(key)
            found_md5 = self.calculate_resource_md5(key)
            if found_md5 != expected_md5:
                return ("File %s referred to by key %s did "
                        "not match its MD5 sum (found: %s, expected %s)") % \
                    (filename, key, found_md5, expected_md5)
        return False

    def _check_refpkg(self):
        error = self.is_invalid()
        if error:
            raise ValueError("Refpkg is invalid: %s" % error)

    @transaction
    def update_metadata(self, key, value):
        """Set *key* in the metadata to *value*.

        Returns the previous value of *key*, or None if the key was
        not previously set.
        """
        old_value = self.contents['metadata'].get(key)
        self.contents['metadata'][key] = value
        self._log('Updated metadata: %s=%s' % (key, value))
        return old_value

    @transaction
    def update_file(self, key, new_path):
        """Insert file *new_path* into the refpkg under *key*.

        The filename of *new_path* will be preserved in the refpkg
        unless it would conflict with a previously existing file, in
        which case a suffix is appended which makes it unique.  The
        previous file, if there was one, is left in the refpkg.  If
        you wish to delete it, see the ``strip`` method.

        The full path to the previous file referred to by *key* is
        returned, or ``None`` if *key* was not previously defined in
        the refpkg.
        """
        if key in self.contents['files']:
            old_path = self.resource_path(key)
        else:
            old_path = None
        self._add_file(key, new_path)
        with open(new_path, 'rb') as f:
            md5_value = md5file(f)
        self.contents['md5'][key] = md5_value
        self._log('Updated file: %s=%s' % (key, new_path))
        if key == 'tree_stats' and old_path:
            warnings.warn('Updating tree_stats, but not phylo_model.',
                          DerivedFileNotUpdatedWarning, stacklevel=2)
        return old_path

    @transaction
    def reroot(self, rppr=None, pretend=False):
        """Reroot the phylogenetic tree.

        This operation calls ``rppr reroot`` to generate the rerooted
        tree, so you must have ``pplacer`` and its auxiliary tools
        ``rppr`` and ``guppy`` installed for it to work.  You can
        specify the path to ``rppr`` by giving it as the *rppr*
        argument.

        If *pretend* is ``True``, the convexification is run, but the
        refpkg is not actually updated.
        """
        with scratch_file(prefix='tree', suffix='.tre') as name:
            # Use a specific path to rppr, otherwise rely on $PATH
            subprocess.check_call([rppr or 'rppr', 'reroot',
                                   '-c', self.path, '-o', name])
            if not(pretend):
                self.update_file('tree', name)
        self._log('Rerooting refpkg')

    def update_phylo_model(self, stats_type, stats_file, frequency_type=None):
        """Parse a stats log and use it to update ``phylo_model``.

        ``pplacer`` expects its input to include the deatils of the
        phylogenetic model used for creating a tree in JSON format
        under the key ``phylo_model``, but no program actually outputs
        that format.

        This function takes a log generated by RAxML or FastTree, parses it,
        and inserts an appropriate JSON file into the refpkg. The first
        parameter must be 'RAxML', 'PhyML' or 'FastTree', depending on which
        program generated the log. It may also be None to attempt to guess
        which program generated the log.

        :param stats_type: Statistics file type. One of 'RAxML', 'FastTree', 'PhyML'
        :param stats_file: path to statistics/log file
        :param frequency_type: For ``stats_type == 'PhyML'``, amino acid
         alignments only: was the alignment inferred with ``model`` or
         ``empirical`` frequencies?

        """

        if frequency_type not in (None, 'model', 'empirical'):
            raise ValueError(
                'Unknown frequency type: "{0}"'.format(frequency_type))
        if frequency_type and stats_type not in (None, 'PhyML'):
            raise ValueError('Frequency type should only be specified for '
                             'PhyML alignments.')

        if stats_type is None:
            with open(stats_file) as fobj:
                for line in fobj:
                    if line.startswith('FastTree'):
                        stats_type = 'FastTree'
                        break
                    elif (line.startswith('This is RAxML') or
                          line.startswith('You are using RAxML')):
                        stats_type = 'RAxML'
                        break
                    elif line.startswith('RAxML-NG'):
                        stats_type = 'RAxML-NG'
                        break
                    elif 'PhyML' in line:
                        stats_type = 'PhyML'
                        break
                else:
                    raise ValueError(
                        "couldn't guess log type for %r" % (stats_file,))

        if stats_type == 'RAxML':
            parser = utils.parse_raxml
        elif stats_type == 'RAxML-NG':
            parser = utils.parse_raxmlng
        elif stats_type == 'FastTree':
            parser = utils.parse_fasttree
        elif stats_type == 'PhyML':
            parser = functools.partial(utils.parse_phyml,
                                       frequency_type=frequency_type)
        else:
            raise ValueError('invalid log type: %r' % (stats_type,))

        with scratch_file(prefix='phylo_model', suffix='.json') as name:
            with open(name, 'w') as phylo_model, open(stats_file) as h:
                json.dump(parser(h), phylo_model, indent=4)
            self.update_file('phylo_model', name)

    def rollback(self):
        """Revert the previous modification to the refpkg.
        """
        # This is slightly complicated because of Python's freakish
        # assignment semantics and because we don't store multiple
        # copies of the log.
        if self.contents['rollback'] is None:
            raise ValueError("No operation to roll back on refpkg")
        future_msg = self.contents['log'][0]
        rolledback_log = self.contents['log'][1:]
        rollforward = copy.deepcopy(self.contents)
        rollforward.pop('rollback')
        self.contents = self.contents['rollback']
        self.contents['log'] = rolledback_log
        self.contents['rollforward'] = [future_msg, rollforward]
        self._sync_to_disk()

    def rollforward(self):
        """Restore a reverted modification to the refpkg.
        """
        if self.contents['rollforward'] is None:
            raise ValueError("No operation to roll forward on refpkg")
        new_log_message = self.contents['rollforward'][0]
        new_contents = self.contents['rollforward'][1]
        new_contents['log'] = [new_log_message] + self.contents.pop('log')
        self.contents['rollforward'] = None
        new_contents['rollback'] = copy.deepcopy(self.contents)
        new_contents['rollback'].pop('rollforward')
        self.contents = new_contents
        self._sync_to_disk()

    def strip(self):
        """Remove rollbacks, rollforwards, and all non-current files.

        When distributing a refpkg, you probably want to distribute as
        small a one as possible.  strip removes everything from the
        refpkg which is not relevant to its current state.
        """
        self._sync_from_disk()
        current_filenames = set(self.contents['files'].values())
        all_filenames = set(os.listdir(self.path))
        to_delete = all_filenames.difference(current_filenames)
        to_delete.discard('CONTENTS.json')
        for f in to_delete:
            self._delete_file(f)
        self.contents['rollback'] = None
        self.contents['rollforward'] = None
        self.contents['log'].insert(
            0, 'Stripped refpkg (removed %d files)' % len(to_delete))
        self._sync_to_disk()

    def start_transaction(self):
        """Begin a transaction to group operations on the refpkg.

        All the operations until the next call to
        ``commit_transaction`` will be recorded as a single operation
        for rollback and rollforward, and recorded with a single line
        in the log.
        """
        if self.current_transaction:
            raise ValueError("There is already a transaction going")
        else:
            initial_state = copy.deepcopy(self.contents)
            self.current_transaction = {'rollback': initial_state,
                                        'log': '(Transaction left no log message)'}

    def commit_transaction(self, log=None):
        """Commit a transaction, with *log* as the log entry."""
        self.current_transaction['rollback'].pop('log')
        self.current_transaction['rollback'].pop('rollforward')
        self.contents['log'].insert(
            0, log and log or self.current_transaction['log'])
        self.contents['rollback'] = self.current_transaction['rollback']
        self.contents['rollforward'] = None  # We can't roll forward anymore
        self.current_transaction = None
        self._sync_to_disk()

    def is_ill_formed(self):
        """Stronger set of checks than is_invalid for Refpkg.

        Checks that FASTA, Stockholm, JSON, and CSV files under known
        keys are all valid as well as calling is_invalid.  Returns
        either False or a string describing the error.
        """
        m = self.is_invalid()
        if m:
            return m

        required_keys = ('aln_fasta', 'aln_sto', 'seq_info', 'tree',
                         'taxonomy', 'phylo_model')
        for k in required_keys:
            if k not in self.contents['files']:
                return "RefPkg has no key " + k

        # aln_fasta, seq_info, tree, and aln_sto must be valid FASTA,
        # CSV, Newick, and Stockholm files, respectively, and describe
        # the same sequences.

        with self.open_resource('aln_fasta') as f:
            firstline = f.readline()
            if firstline.startswith('>'):
                f.seek(0)
            else:
                return 'aln_fasta file is not valid FASTA.'
            fasta_names = {seq.id for seq in fastalite(f)}

        with self.open_resource('seq_info') as f:
            lines = list(csv.reader(f))
            headers = set(lines[0])

            # Check required headers
            for req_header in 'seqname', 'tax_id':
                if req_header not in headers:
                    return "seq_info is missing {0}".format(req_header)
            lengths = {len(line) for line in lines}
            if len(lengths) > 1:
                return "some lines in seq_info differ in field cout"
            csv_names = {line[0] for line in lines[1:]}

        with self.open_resource('aln_sto') as f:
            try:
                sto_names = set(utils.parse_stockholm(f))
            except ValueError:
                return 'aln_sto file is not valid Stockholm.'

        try:
            tree = dendropy.Tree.get(
                path=self.resource_path('tree'),
                schema='newick',
                case_sensitive_taxon_labels=True,
                taxon_namespace=dendropy.TaxonNamespace(is_case_sensitive=True),
                preserve_underscores=True)
            tree_names = set(tree.taxon_namespace.labels())
        except Exception:
            return 'tree file is not valid Newick.'

        d = fasta_names.symmetric_difference(sto_names)
        if len(d) != 0:
            return "Names in aln_fasta did not match aln_sto.  Mismatches: " + \
                ', '.join([str(x) for x in d])
        d = fasta_names.symmetric_difference(csv_names)
        if len(d) != 0:
            return "Names in aln_fasta did not match seq_info.  Mismatches: " + \
                ', '.join([str(x) for x in d])
        d = fasta_names.symmetric_difference(tree_names)
        if len(d) != 0:
            return "Names in aln_fasta did not match nodes in tree.  Mismatches: " + \
                ', '.join([str(x) for x in d])

        # Next make sure that taxonomy is valid CSV, phylo_model is valid JSON
        with self.open_resource('taxonomy') as f:
            lines = list(csv.reader(f))
            lengths = {len(line) for line in lines}
            if len(lengths) > 1:
                return ("Taxonomy is invalid: not all lines had "
                        "the same number of fields.")
            # I don't try to check if the taxids match up to those
            # mentioned in aln_fasta, since that would make taxtastic
            # depend on RefsetInternalFasta in romperroom.
        with self.open_resource('phylo_model') as f:
            try:
                json.load(f)
            except ValueError:
                return "phylo_model is not valid JSON."

        return False

    def load_db(self):
        """Load the taxonomy into a sqlite3 database.

        This will set ``self.db`` to a sqlite3 database which contains all of
        the taxonomic information in the reference package.
        """

        db = taxdb.Taxdb()
        db.create_tables()
        reader = csv.DictReader(self.open_resource('taxonomy', 'rU'))
        db.insert_from_taxtable(lambda: reader._fieldnames, reader)

        curs = db.cursor()
        reader = csv.DictReader(self.open_resource('seq_info', 'rU'))
        curs.executemany("INSERT INTO sequences VALUES (?, ?)",
                         ((row['seqname'], row['tax_id']) for row in reader))

        db.commit()
        self.db = db

    def most_recent_common_ancestor(self, *ts):
        """Find the MRCA of some tax_ids.

        Returns the MRCA of the specified tax_ids, or raises ``NoAncestor`` if
        no ancestor of the specified tax_ids could be found.
        """
        if len(ts) > 200:
            res = self._large_mrca(ts)
        else:
            res = self._small_mrca(ts)

        if res:
            (res,), = res
        else:
            raise NoAncestor()
        return res

    def _large_mrca(self, ts):
        """Find the MRCA using a temporary table."""
        cursor = self.db.cursor()

        cursor.execute("""
            DROP TABLE IF EXISTS _mrca_temp
        """)

        cursor.execute("""
            CREATE TEMPORARY TABLE _mrca_temp(
                child TEXT PRIMARY KEY REFERENCES taxa (tax_id) NOT NULL
            )
        """)

        cursor.executemany("""
            INSERT INTO _mrca_temp
            VALUES (?)
        """, ((tid,) for tid in ts))
        cursor.execute("""
            SELECT parent
            FROM   _mrca_temp
                   JOIN parents USING (child)
                   JOIN taxa
                     ON parent = taxa.tax_id
                   JOIN ranks USING (rank)
            GROUP  BY parent
            HAVING COUNT(*) = ?
            ORDER  BY rank_order DESC
            LIMIT  1
        """, (len(ts),))

        return cursor.fetchall()

    def _small_mrca(self, ts):
        """Find a MRCA using query parameters.

        This only supports a limited number of tax_ids; ``_large_mrca`` will
        support an arbitrary number.
        """
        cursor = self.db.cursor()
        qmarks = ', '.join('?' * len(ts))
        cursor.execute("""
            SELECT parent
            FROM   parents
                   JOIN taxa
                     ON parent = taxa.tax_id
                   JOIN ranks USING (rank)
            WHERE  child IN (%s)
            GROUP  BY parent
            HAVING COUNT(*) = ?
            ORDER  BY rank_order DESC
            LIMIT  1
        """ % qmarks, ts + (len(ts),))

        return cursor.fetchall()

    def file_abspath(self, resource):
        """Deprecated alias for *resource_path*."""
        warnings.warn(
            "file_abspath is deprecated; use resource_path instead",
            DeprecationWarning, stacklevel=2)
        return self.resource_path(resource)

    def file_name(self, resource):
        """Deprecated alias for *resource_name*."""
        warnings.warn(
            "file_name is deprecated; use resource_name instead",
            DeprecationWarning, stacklevel=2)
        return self.resource_name(resource)

    def file_md5(self, resource):
        """Deprecated alias for *resource_md5*."""
        warnings.warn(
            "file_md5 is deprecated; use resource_md5 instead",
            DeprecationWarning, stacklevel=2)
        return self.resource_md5(resource)
