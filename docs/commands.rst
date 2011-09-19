Detailed ``taxit`` documentation
================================

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
