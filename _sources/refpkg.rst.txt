Creating and manipulating refpkgs
=================================

Quickstart
----------

The refpkgs is a container format for keeping all the miscellaneous files required to run pplacer_ in one place.  Roughly, it is a directory containing files and JSON file describing the other contents.  Refpkgs can be created and manipulated either from the commandline with subcommands of ``taxit``, or from Python using an API exposed by the module ``taxtastic.refpkg``.

To create an empty refpkg from Python, you would write::

    from taxtastic.refpkg import Refpkg

    r = Refpkg('/path/to/new/refpkg')

If there had already been a refpkg at ``/path/to/new/refpkg``, ``r`` would contain a reference to the existing refpkg afterwards.  If it doesn't exist, it is created.  For historical reasons, the command line interface requires that you specify a locus to be recorded in the refpkg's metadata::

    taxit create -l locus_name -P /path/to/new/refpkg

``taxit create`` takes many other arguments to add particular files that ``pplacer`` expects, but we won't go into them here.  To add a file to a refpkg, we specify a key which it will be added under (refpkgs act as key-value stores) and the file to add.  With the API, this is::

    r.update_file('key', '/path/to/file/to/add')

or from the command line

    taxit update /path/to/refpkg key=/path/to/file/to/add

Either way, ``/path/to/file/to/add`` is copied into the refpkg (and renamed if necessary so as not to collide with files already in the refpkg), and an entry added to the JSON file to assign ``key`` to the file.

Refpkgs also store metadata, in the form of a key-value store of arbitrary strings.  The metadata keys occupy a separate namespace from the file keys, so you can have both a file and a string assigned to the key ``boris`` at the same time.  Setting metadata with the API uses the method ``update_metadata``::

    r.update_metadata('key', 'value to add')

From the command line, use ``taxit update``, but with the ``--metadata`` option::

    taxit update --metadata /path/to/refpkg "key=value to add"

Refpkgs have an undo/redo mechanism.  If you update something incorrectly, you can call the ``rollback`` method of the API or the ``rollback`` subcommand of ``taxit`` to restore the refpkg to its previous state.  An improperly rolled back operation can be rolled forward again with ``rollforward``.  So running the following in the API leaves the refpkg unchanged at the end::

    r.rollback()
    r.rollforward()

Similarly with ``taxit``::

    taxit rollback /path/to/refpkg
    taxit rollforward /path/to/refpkg

Before distributing a refpkg, you may want to strip out rollback/rollforward information and any unused files.  The ``strip`` method removes everything not necessary to the current state of the refpkg::

    r.strip()

or from the command line::

    taxit strip /path/to/refpkg

There is a method and a command, both called ``check``, which check if a refpkg is usable as an input to ``pplacer``.  Running ``r.check()`` from Python or ``taxit check /path/to/refpkg`` from the command line will fail, saying that there is no such key ``aln_fasta`` in the refpkg.

The Refpkg format
-----------------

A refpkg is a crude key-value store for files with some basic integrity checking, machinery for undo and redo of operations, and a little bit of metadata storage.  The minimal refpkg consists of a directory containing a file names ``CONTENTS.json``.  The JSON file must contain the following keys:

``files``
  A list of the files this refpkg is currently tracking in its directory.   The value must be a JSON object assigning keys to filenames in the directory, e.g. ``{"taxonomy": "taxtable.csv"}``.
``md5``
  The value must be a JSON object with the same keys are ``files``, but where the values are the MD5 sums of the files.  These fields are used to ensure the integrity of the files the refpkg is tracking.
``metadata``
  The value must be a JSON object containing keys referring to strings, e.g., ``{"author": "Boris the mad baboon", "create_date": "2011-08-18 14:50:39"}``.  This is where any data describing the refpkg as a whole should be.
``log``
  The value must be a list of strings, e.g., ``["Created package.", "Oh god, get it off me!"]``.  The various refpkg operations each append an entry to the log, so it records the history of what has been done to this refpkg.
``rollback``
  Either ``null`` or the JSON object which was previously the top level object of ``CONTENTS.json`` before the last operation performed on this refpkg.  It is used to undo operations on a refpkg.
``rollforward``
  When an operation is rolled back, the state before the rollback is preserved in ``rollforward`` so the undo can be redone.  ``rollforward`` is either ``null`` or a list of two entries, the first a string giving the log entry associated with the rolled back operation, the second the JSON object describing the contents before the rollback.

Any program only wanting to read refpkgs only needs to worry about the keys ``files``, ``md5``, and ``metadata``.  Any file read from the refpkg should have its MD5 sum checked against the refpkg's stored value.

The refpkg format was designed to store multiple alignments and trees with optional taxonomic information for use by ``pplacer``, so certain fields are expected.

``taxonomy``
  A CSV file giving a taxonomy that contains all the sequences included in the refpkg.  This is typically created with taxtastic's ``taxit taxtable`` command.

``profile``
  The profile used by ``cmalign`` or ``hmmalign`` when aligning the sequences for ``pplacer``.

``tree``
  A phylogenetic tree with the sequences in this refpkg as its nodes, stored in Newick format.

``tree_stats``
  The output of the phylogenetic inference program describing its run when it assembled the phylogenetic tree for this refpkg.

``phylo_model``
  A JSON file describing the phylogenetic model used for tree construction, usually parsed from the information in ``tree_stats``.

``aln_fasta``
  A FASTA file containing all the sequences included in this refpkg.

``seq_info``
  A CSV file giving basic information on all the sequences included in the refpkg.  It should begin with one line giving the field names::

      "seqname","accession","tax_id","species_name","is_type"

  ``seqname`` should match the ID of the sequence in the FASTA file.  ``accession`` is a database reference for the sequence, which can be the same as ``seqname``.  In our work, ``seqname`` is an RDP accession number and ``accession`` is the NCBI accession number corresponding to that RDP entry.  ``tax_id`` is the entry in the taxonomy this sequence is mapped to, and ``species_name`` is the name associated with that entry.  ``is_type`` indicates whether this sequence is from a typestrain or not (again, this is particular to our work).

``aln_sto``
  The same sequences as in ``aln_fasta``, but written in Stockholm format.

The files referred to by ``aln_fasta``, ``seq_info``, ``aln_sto``, and ``tree`` should all have the same list of sequences.  This isn't strictly enforced, but you can check that it is so with the ``taxit check`` command or the ``is_ill_formed`` method of the refpkg API.

The undo/redo system is implemented as a purely functional data structure known as a zipper, used as a replacement for arrays with a pointer into them as a cursor in languages where data is immutable.  A zipper consists of a current entry, a ordered list of previous entries, and an ordered list of subsequent entries.  Moving the cursor one place to the right is equivalent to pushing the current entry onto the head of the list of previous entries, and popping the head of the list of subsequent entries and making that the new current entry.  Moving the cursor one place to the left is exactly the opposite.

The top level object of ``CONTENTS.json`` plays the role of the current entry, and its fields ``rollback`` and ``rollforward`` are the heads of the lists of previous and subsequent entries.  The ``rollback`` field of the object in the ``rollback`` field is the second element of the list of previous entries, etc.  Thus undoing an operation consists of putting the current toplevel JSON object in the ``rollforward`` fields of the object in the ``rollback`` field, and making that object the new toplevel JSON object (with some book keeping details to keep everything consistent).

As a result of this, there may be files besides those referenced in the ``files`` key of the JSON object in the refpkg.  They may be referenced by other entries in the zipper.  There is no attempt to intelligently garbage collect orphaned files.  They are only deleted when the refpkg's ``strip`` method is called, which removes all undo/redo information as well.



.. _pplacer: http://matsen.fhcrc.org/pplacer
