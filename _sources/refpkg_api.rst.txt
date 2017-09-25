The Python Refpkg API
=====================

Taxtastic provides a Python API for creating and manipulating refpkgs wihtout having to recourse to the command line.  It is entirely based around the ``Refpkg`` class in ``taxtastic.refpkg``.  Thus all scripts dealing with refpkgs will include code that looks something like::

    from taxtastic.refpkg import Refpkg

    r = Refpkg('/path/to/refpkg')

The constructor takes no other arguments.  If the refpkg at ``/path/to/refpkg`` already exists, ``r`` is attached to it.  If ``/path/to/refpkg`` does not exist, a new, empty refpkg is created.  If there is already something at ``/path/to/refpkg``, but it is not a refpkg, ``Refpkg`` throws a ``ValueError``.

The rest of the API is implemented as methods of the ``Refpkg`` class.

**Note**: This library is *not* thread safe!

Accessing refpkg contents
-------------------------

Refpkgs are primarily key-value stores for files and metadata.  To find all the keys referring to files or to metadata in a given refpkg, use the methods

.. automethod:: taxtastic.refpkg.Refpkg.file_keys

.. automethod:: taxtastic.refpkg.Refpkg.metadata_keys

For example::

    r = Refpkg('/path/to/new/refpkg')
    r.metadata_keys()
     --> ['create_date', 'format_version']
    r.update_file('file_key', '/path/to/some/file')
    r.file_keys()
     --> ['file_key']

Call the ``metadata`` method to retrieve the value of a particular metadata field.

.. automethod:: taxtastic.refpkg.Refpkg.metadata

For instance, on the same refpkg as above::

    r.metadata('format_version')
     --> '1.1'

There are several methods for working with the resources in a reference package.

.. automethod:: taxtastic.refpkg.Refpkg.open_resource

.. automethod:: taxtastic.refpkg.Refpkg.resource_name

.. automethod:: taxtastic.refpkg.Refpkg.resource_md5

.. automethod:: taxtastic.refpkg.Refpkg.resource_path

Checking refpkg integrity
-------------------------

There are two methods for checking a refpkg.  The ``is_invalid`` method enforces only that the refpkg is sane: there is a ``CONTENTS.json`` file, the MD5 sums listed match the actual files, and other such basics.  The ``is_ill_formed`` method is much stronger.  It enforces the necessary structure of a refpkg to be fed to pplacer_.

.. automethod:: taxtastic.refpkg.Refpkg.is_invalid

.. automethod:: taxtastic.refpkg.Refpkg.is_ill_formed

Updating and modifying refpkgs
------------------------------

All updates are made via two methods: ``update_metadata`` and ``update_file``.  Each takes a key and a new value to set at that key, and returns the path to the previous value of the key, or ``None`` if the key was not previously defined in the refpkg.

.. automethod:: taxtastic.refpkg.Refpkg.update_metadata

.. automethod:: taxtastic.refpkg.Refpkg.update_file


Refpkg history, undo, and redo
------------------------------

Each operation performed on a refpkg leaves an entry in a log stored in ``CONTENTS.json``.  You can access this log by calling the ``log`` method.

.. automethod:: taxtastic.refpkg.Refpkg.log

Each operation call also be undone (and redone once undone).  The undo stack is arbitrarily deep so all operations back to the previous call to ``strip`` (see below) can be undone.  To undo an operation, call ``rollback``.  To redo it afterwards, call ``rollforward``.  The logs will similarly be updated to stay in sync with the rollback and rollforward of operations.  Note that when you call another operation, all the redo information before that point is removed.  You cannot undo an operation, perform another operation, then redo the first operation.

.. automethod:: taxtastic.refpkg.Refpkg.rollback

.. automethod:: taxtastic.refpkg.Refpkg.rollforward

After performing a lot of operations on a refpkg, there will often be a long undo history, and files no longer referred to in the refpkg's current state.  To remove everything not relevant to the refpkg's current state other than the log, call the ``strip`` method.

.. automethod:: taxtastic.refpkg.Refpkg.strip

You can force a series of operations to be recorded as a single operation for rollback and have a single log entry by calling ``start_transaction`` before them, and ``commit_transaction`` with the log entry to record when they are done.

.. automethod:: taxtastic.refpkg.Refpkg.start_transaction

.. automethod:: taxtastic.refpkg.Refpkg.commit_transaction

For example, ::

    from taxtastic.refpkg import Refpkg

    r = Refpkg('/path/to/refpkg')

    r.start_transaction()
    r.update_metadata('author', 'Boris the mad baboon')
    r.update_file('boris_signature', '/path/to/some/file')
    r.commit_transaction("Left Boris's mark!")

would result in a single operation that could be rolled back as one, and leaves the log entry ``"Left Boris's mark!"``.

``pplacer`` specific commands
-----------------------------

Finally, the API has three commands which are specific to creating inputs for pplacer_.  One of these is ``check``, which was described above.  The other two are:

.. automethod:: taxtastic.refpkg.Refpkg.reroot

.. automethod:: taxtastic.refpkg.Refpkg.update_phylo_model

.. _pplacer: http://matsen.fhcrc.org/pplacer
