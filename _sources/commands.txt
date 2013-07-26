Detailed ``taxit`` documentation
================================

This section gives the detailed documentation on ``taxit``'s subcommands, organized alphabetically.


add_nodes
---------

``taxit add_nodes [...] -d database_file -N nodes.csv``

Add additional nodes, specified in ``nodes.csv`` to the taxonomy in ``database_file``.  ``nodes.csv`` should be a CSV file with its first line naming fields or an Excel 97 (``.xls``) file.  It *must* specify the columns

``tax_id``
  The tax_id for this new node in the taxonomy, which must not conflict with an existing tax_id.  NCBI's tax_ids are all integers, so it works well to choose an alphabetic prefix for the tax_ids for all new nodes (e.g., name them AB1, AB2, AB3, etc.).
``rank``
  The name of the rank at which this node falls in the taxonomy (e.g., ``species`` or ``genus``).  The ranks understood by the NCBI taxonomy are ``root``, ``superkingdom``, ``kingdom``, ``subkingdom``, ``superphylum``, ``phylum``, ``subphylum``, ``superclass``, ``class``, ``subclass``, ``infraclass``, ``superorder``, ``order``, ``suborder``, ``parvorder``, ``infraorder``, ``superfamily``, ``family``, ``subfamily``, ``tribe``, ``subtribe``, ``genus``, ``subgenus``, ``species``, ``group``, ``species``, ``subgroup``, ``species``, ``subspecies``, ``varietas``, and ``forma``.
``parent_id``
  The tax_id which will be set as the parent of this node in the taxonomy.
``tax_name``
  The name of this taxon.

Any combination of the following columns *may* be specified:

``source_name``
  A string describing the origin of these taxa so that it is easy to find them in the database.
``source_id``
  An integer specifying the source of these new nodes, used internally in the database to match to ``source_name``.
``children``
  A semicolon separated list of tax_ids which should be detached from their current parents and attached to this node as its children.

Examples::

    # Add the taxa specified in new_taxa.csv to taxonomy.db
    taxit add_nodes -d taxonomy.db -N new_taxa.csv

    # Add the same taxa, but force the source name to be hetty_lab
    taxit add_nodes -d taxonomy.db -N new_taxa.csv -S hetty_lab

Arguments:

``-h``
  Print help and exit.
``-d``, ``--database-file``
  Add the new new nodes to the specified database, which should be one created by the ``taxit new_database`` command.
``-N``, ``--new-nodes``
  Pull the new nodes from this file.  The format is specified above.
``-S``, ``--source-name``
  Force all the nodes added from this command to have a particular source name.

check
-----

``taxit check /path/to/refpkg``

Check whether ``/path/to/refpkg`` is a valid input for ``pplacer``, that is, does it have a FASTA file of the reference sequences, a Stockholm file of their multiple alignment, a Newick formatted tree build from the aligned sequences, and all the necessary auxiliary information.


create
------

``taxit create [...] -P refpkg -l "Locus description"``

Create a new refpkg at the location specified by the argument to ``-P`` with locus name ``-l``.  All other fields are used to specify initial metadata and files to add to the refpkg.  If there is already a refpkg at ``refpkg``, this command will fail unless you specify ``-c`` or ``--clobber``.

**General options**

``-c``, ``--clobber``
  If there is an existing refpkg at the given path, delete it and create a new one.
``-P``, ``--package-name``
  Specify the refpkg to create (required).

**Metadata options**

=============================     ===================
          Option                  Metadata key
=============================     ===================
``-a``, ``--author``              ``author``
``-d``, ``--description``         ``description``
``-l``, ``--locus``               ``locus``
``-r``, ``--package-version``     ``package_version``
=============================     ===================


**File options**

=========================     ==============     =====================================
          Option                File key         Description
=========================     ==============     =====================================
``-f``, ``--aln-fasta``       ``aln_fasta``      Reference sequences in FASTA format
``-i``, ``--seq-info``        ``seq_info``       CSV describing aligned sequences
``-m``, ``--mask``            ``mask``           Text file containing sequence mask
``-p``, ``--profile``         ``profile``        Multiple alignment profile
``-R``, ``--readme``          ``readme``         A README file for the refpkg
``-s``, ``--tree-stats``      ``tree_stats``     Typically written by the tree builder
``-S``, ``--aln-sto``         ``aln_sto``        Stockholm file of reference sequences
``-t``, ``--tree-file``       ``tree``           Phylogenetic tree in Newick format
``-T``, ``--taxonomy``        ``taxonomy``       CSV file specifying taxonomy
=========================     ==============     =====================================

Examples::

    # Create a minimal refpkg
    taxit create -P my_refpkg -l "Some locus name"

    # Create a refpkg with lots of files in it
    taxit create -P another_refpkg -l "Another locus" \
        --author "Boris the mad baboon" --package-version 0.3.1 \
        --aln-fasta seqs.fasta --aln-sto seqs.sto \
        --tree-file seqs.newick --seq-info seqs.csv \
        --profile cmalign.profile --tree-stats RAxML.info \
        --taxonomy taxtable.csv


findcompany
-----------

``taxit findcompany [-c] taxonomy.db [-i taxids.txt] [taxid ...] [-o output.txt]``

A command meant to follow ``lonelynodes`` (below). Given a list of tax_ids produced by ``taxit lonelynodes``, produces another list of species tax_ids that can be added to the taxtable that would render those tax_ids no longer lonely.

Examples::

    taxit findcompany taxonomy.db -i taxids.txt -o newtaxids.txt
    taxit findcompany taxonomy.db 31661 5213 564

Arguments::

``-c``
    Produce only one output tax_id per input tax_id, whether or not the output species would themselves be lonely.

``-o``
    Write new taxids to the specified file. Otherwise they are written to ``stdout``

``-i``
    Read taxids from the specified file in addition to any given as command line arguments.


lonelynodes
-----------

``taxit lonelynodes target [-o output.txt]``

Find nodes in ``target`` (which can be a CSV file extracted by ``taxit taxtable`` or a RefPkg containing such a file) which are lonely -- that is, whose parents have only one child. Print them, one per line, to ``stdout`` or to the file specified by the ``-o`` option.

Examples::

    # Find lonely nodes in RefPkg mypkg-0.1.refpkg
    taxit lonelynodes mypkg-0.1.refpkg

Arguments::

``-h``
  Print help and exit

``-o``
  Write resulting tax_ids to a specified filename instead of ``stdout``.

``-v``
  Run verbosely.

new_database
------------

``taxit new_database [...] -d database_file``

Download the current version of the NCBI taxonomy and load it into ``database_file`` as an SQLite3 database.  If ``database_file`` already exists, it will fail and leave it untouched unless you specify ``-x`` or ``--clobber``.  The NCBI taxonomy will be downloaded into the same directory as ``database_file`` will be created in unless you specify ``-p`` or ``--download-dir``.

Examples::

    # Download the NCBI taxonomy and create taxonomy.db if it does not exist
    taxit new_database -d taxonomy.db
    
    # Force the creation of taxonomy.db in the parent directory, putting
    # the downloaded NCBI data in /tmp/ncbi.
    taxit new_database -d ../taxonomy.db -x -p /tmp/ncbi

Arguments:

``-h``
  Print help and exit.

``-d``, ``--database-file``
  Specify the file in which the NCBI taxonomy should be written.

``--x``, ``--clobber``
  Replace ``database_file`` if it already exists.

``-p``, ``--download-dir``
  Download the NCBI taxonomy into the specified path.  If not specified, the taxonomy will be downloaded into the same directory where the final database will be created.

reroot
------

``taxit reroot refpkg``

Calls ``rppr reroot`` to generate a rerooted tree from the tree in ``refpkg`` and writes it back to the refpkg.  The refpkg ``refpkg`` must contain the necessary inputs for ``pplacer`` for this to work.

Examples::

    # Reroot the tree in my_refpkg
    taxit reroot my_refpkg

    # Try running reroot without modifying the refpkg, using a particular 
    # version of rppr
    taxit reroot --rppr ~/local/bin/rppr -p my_refpkg

Arguments:

``--rppr``
  Specify the path to the ``rppr`` executable to use.
``-p``, ``--pretend``
  Calculate the rerooted tree, but don't actually change the tree file in the refpkg.


rollback
--------

``taxit rollback [-n N] refpkg``

Rollback ``N`` operations on ``refpkg`` (default to 1 operation if ``-n`` is omitted).  This is equivalent to calling the ``rollback()`` method of ``taxtastic.refpkg.Refpkg``.  If there are not at least ``N`` operations that can be rolled back, an error is returned and no changes are made to the refpkg.

Examples::

    # Update the author on my_refpkg, then roll back the change so
    # that it is in the same state it was.
    taxit update --metadata 'author=Boris the mad baboon'
    taxit rollback my_refpkg

    # Roll back the last 3 operations on my_refpkg
    taxit rollback -n 3 my_refpkg

Arguments:

``-n``
  Give an integer specifying the number of operations to roll back (default: 1)


rollforward
-----------

``taxit rollforward [-n N] refpkg``

Restore the last ``N`` rolled back operations on ``refpkg``, or the last operation if ``-n`` is omitted.  If there are not at least ``N`` operations that can be rolled forward on this refpkg, then an error is returned and no changes are made to the refpkg.

Note that operations can only be rolled forward immediately after being rolled back.  If any operation besides a rollback occurs, all roll forward information is removed.

Examples::

    # Roll back the last operation on my_refpkg, then restore it.
    taxit rollback my_refpkg
    taxit rollforward my_refpkg

    # Roll forward the last 3 rollbacks on my_refpkg
    taxit rollforward -n 3 my_refpkg

Arguments:

``-n``
  Give an integer specifying the number of operations to roll back (default: 1)



strip
-----

``taxit strip refpkg``

Delete everything in the refpkg not relevant to the current state: all files not referred to in the current state and all rollback and rollforward information.  The log is preserved, with a new entry entered indicating that ``refpkg`` was stripped.

Examples::

    taxit update my_refpkg hilda=file1

    # After this, file1 is still in the refpkg, but not referred to
    # except by the rollback information.
    taxit update my_refpkg hilda=file2

    # strip deletes file1, and the rollback and rollforward information
    taxit strip my_refpkg


taxids
------

``taxit taxids -d ncbi_taxonomy.db [-f file_of_names.txt|-n name1,name2,...] [-o taxids.txt]``

Convert a list of taxonomic names into a list of tax_ids.  ``ncbi_taxonomy.db`` must be a database created by ``taxit new_database``, containing a taxonomy.  The names to convert can be specified in a text file with one name per line (the ``-f`` or ``--name-file`` options) or on the command line as a comma delimited list (the ``-n`` of ``--name`` options).

Examples::

    # Look up two species and print their tax_ids to stdout, one per line
    taxit taxids -d ncbi_database.db -n "Lactobacillus crispatus,Lactobacillus helveticus"

    # Read the species from some_names.txt and write their tax_ids to some_taxids.txt
    taxit taxids -d ncbi_database.db -f some_names.txt -o some_taxids.txt

Arguments:

``-d``, ``--database-file``
  The taxonomy file to look up names in.  Must be an SQLite3 database as created by ``taxit new_database``. (Required)
``-f``, ``--name-file``
  Read names, one per line, from this file to convert to tax_ids.
``-n``, ``--name``
  Specify a comma separated list of names to look up in the taxonomy and convert to tax_ids.
``-o``, ``--out-file``
  Write the tax_ids looked up in the taxonomy to this file.  (Default: stdout)


taxtable
--------

``taxit taxtable [...] -d database_file [-n taxa_names.txt] [-t tax_ids] [-o output.csv]``

Write a CSV file containing the minimal subset of the taxonomy in ``database_file`` which encompasses all the taxa specified in ``taxa_names.txt`` and ``tax_ids`` and all nodes connecting them to the root of the taxonomy.  By default the CSV is written to ``stdout``, unless redirectored with ``-o`` or ``--out-file``.

``taxa_names.txt`` should be a text file specifying names of taxa.  Python style comments are ignored as are empty lines.  Names may be separated by commas, semicolons, and arbitrary amounts of whitespace on both sides of those separators, but the whitespace within them must be exact (e.g., ``Lactobacillus crispati`` must have exactly one space between the two words to match the entry in the taxonomy).

``tax_ids`` is either a comma or semicolon delimited list of tax_ids (e.g., ``4522,2213;44;221``) or the name of a text file containing tax_ids.  The text file also allows Python style comments, and any non-comment text separated by and combination of spaces, commas, and semicolons is considered a tax_id.

tax_ids and taxa names can overlap, nor does anything have to be unique within either file.  The nodes will only be written once in the CSV output no matter how many times a particular taxon is mentioned.

Examples::

    # Extract tax_ids 47770 and 33945 and all nodes connecting them to the root.
    taxit taxtable -d taxonomy.db -t 47770,33945

    # The same as above, but write the output to subtax.csv instead of stdout
    taxit taxtable -d taxonomy.db -t 47770,33945 -o subtax.csv

    # Extract the same tax_ids, plus the taxa specifies in taxnames.txt
    taxit taxtable -d taxonomy.db -t 47770,33945 -n taxnames.txt -o taxonomy_from_both.csv

Arguments:

``-h``
  Print help and exit.
``-d``, ``--database-file``
  Use the specified database as the taxonomy to subset.  The database should be one created by ``taxit new_database``.
``-n``, ``--tax-names``
  Include these taxa names and all nodes connecting them to the root of the taxonomy in the output.
``-t``, ``--tax-ids``
  Include these tax_ids and all nodes connecting them to the root of the taxonomy in the output.  The argument can be either a filename or a list of tax_ids separated by commas or semicolons.
``-o``, ``--out-file``
  Write the output to the given filename instead of stdout.


update
------

``taxit update refpkg [--metadata] "key=some value" ...``

Update ``refpkg`` to set ``key`` to ``some value``.  If ``--metadata`` is specified, the update is done to the metadata.  Otherwise ``some value`` is treated as the path to a file, and that file is updated in ``refpkg``.  An arbitrary of "key=value" pairs can be specified on the command line.  If the same key is specified twice, the later occurrence dominates.

All updates specified to an instance of this command are run as a single operation, and will all be undone by a single rollback.

Examples::

    # Set the author in my_refpkg
    taxit update my_refpkg --metadata "author=Boris the mad baboon"

    # Set the author and version at once
    taxit update my_refpkg --metadata "author=Bill" "package_version=1.7.2"

    # Insert a file into the refpkg
    taxit update my_refpkg "aln_fasta=/path/to/a/file.fasta"

Arguments:

``--metadata``
  Treat all the updates as changes to metadata, not files.
