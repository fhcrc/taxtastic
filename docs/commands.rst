Detailed ``taxit`` documentation
================================

This section gives the detailed documentation on ``taxit``'s subcommands, organized alphabetically.


add_nodes
---------

.. literalinclude:: _helptext/add_nodes.txt

Add nodes or names to the taxonomy in the specified database.
``new_nodes`` should be a yaml format file containing one or more
records, each of which specifies a new node or name.

For a new **node** the following are required:

``type``
  The value must be "node"
``tax_id``
  The tax_id for this new node in the taxonomy, which must not
  conflict with an existing tax_id.  NCBI's tax_ids are all integers,
  so it works well to choose an alphabetic prefix for the tax_ids for
  all new nodes (e.g., name them AB1, AB2, AB3, etc.).
``rank``
  The name of the rank at which this node falls in the
  taxonomy. Choose from among the ranks specified in table ``ranks``.
``parent_id``
  The tax_id which will be set as the parent of this node in the taxonomy.
``names``
  One or more names to associate with the node. Minimally, must define
  a single taxonomic name. See description of a name record below.

Required unless a default source name is specified using the
``--source-name`` option:

``source_name``
  A string describing the origin of these taxa so that it is easy to
  find them in the database.

Any combination of the following columns *may* be specified:

``children``
  A list of tax_ids which should be detached from their current
  parents and attached to this node as its children.

A minimal example of a record specifying a node (assuming
``--source-name`` is provided on the command line)::

  ---
  type: node
  tax_id: "newid"
  parent_id: "1279"
  rank: species_group
  names:
    - tax_name: between genus and species

A record providing ``source_name``, multiple taxonomic names, plus
child nodes::

  ---
  type: node
  tax_id: "newid"
  parent_id: "1279"
  rank: species_group
  names:
    - tax_name: between genus and species
      is_primary: true
    - tax_name: another name
  source_name: someplace
  children:
    - "1280" # Staphylococcus aureus
    - "1281" # Staphylococcus carnosus

A record specifying **names** to be added to existing nodes has the
following required fields:

``type``
  The value must be "name"
``tax_id``
  The tax_id to add names to.
``names``
  A list of taxonomic names. If a single name is provided, requires
  only ``tax_name``; if more than one, the primary name must be
  indicated (see examples below).

A minimal example (again, assuming ``source_name`` is defined from the command line)::

  ---
  type: name
  tax_id: bar
  names:
    - tax_name: a new name for bar

If there are multiple names::

  ---
  type: name
  tax_id: bar
  names:
    - tax_name: a new name for bar
      is_primary: true
    - tax_name: another name

Multiple records are delimited by ``---`` and may contain any
combination of names and nodes::

  ---
  type: node
  tax_id: "newid"
  parent_id: "1279"
  rank: species_group
  names:
    - tax_name: between genus and species
  ---
  type: name
  tax_id: bar
  names:
    - tax_name: a new name for bar

Note that the nodes and names are added to the database in the order
specified; be sure to add parent nodes before children.

add_to_taxtable
---------------

.. literalinclude:: _helptext/add_to_taxtable.txt

check
-----

.. literalinclude:: _helptext/check.txt

composition
-----------

.. literalinclude:: _helptext/composition.txt

create
------

.. literalinclude:: _helptext/create.txt

**Input files**

Input files are identified in the refpkg using the following labels
(see, for example ``taxit rp``):

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

.. literalinclude:: _helptext/findcompany.txt

Examples::

  taxit findcompany taxonomy.db -i taxids.txt -o newtaxids.txt
  taxit findcompany taxonomy.db 31661 5213 564

info
----

.. literalinclude:: _helptext/info.txt


lineage_table
-------------

.. literalinclude:: _helptext/lineage_table.txt

Examples::

  taxit taxtable taxonomy.db -i seq_info.csv -o taxtable.csv
  taxit lineage_table taxtable.csv seq_info.csv \
      --csv-table taxonomy.csv --taxonomy-table taxonomy.txt

``taxonomy.txt`` looks like this::

  s1	"pk__Bacteria";"ph__Firmicutes";"cl__Bacilli";"or__Bacillales";"fa__Staphylococcaceae";"ge__Staphylococcus";"sp__Staphylococcus aureus"
  s2	"pk__Bacteria";"ph__Firmicutes";"cl__Bacilli";"or__Bacillales";"fa__Staphylococcaceae";"ge__Staphylococcus";"sp__Staphylococcus equorum"
  s3	"pk__Bacteria";"ph__Firmicutes";"cl__Bacilli";"or__Bacillales";"fa__Staphylococcaceae";"ge__Staphylococcus";"sp__Staphylococcus equorum"
  s4	"pk__Bacteria";"ph__Firmicutes";"cl__Bacilli";"or__Bacillales";"fa__Staphylococcaceae";"ge__Staphylococcus";"sp__unclassified"


lonelynodes
-----------

.. literalinclude:: _helptext/lonelynodes.txt

Examples::

    # Find lonely nodes in RefPkg mypkg-0.1.refpkg
    taxit lonelynodes mypkg-0.1.refpkg

new_database
------------

.. literalinclude:: _helptext/new_database.txt

Examples:

    Download the NCBI taxonomy and create taxonomy.db if it does not exist::

      taxit new_database taxonomy.db

    Force the creation of taxonomy.db in the parent directory, putting
    the downloaded NCBI data in /tmp/ncbi::

      taxit new_database ../taxonomy.db -p /tmp/ncbi

refpkg_intersection
-------------------

.. literalinclude:: _helptext/refpkg_intersection.txt

reroot
------

.. literalinclude:: _helptext/reroot.txt

Examples:

Reroot the tree in my_refpkg::

  taxit reroot my_refpkg

Try running reroot without modifying the refpkg, using a particular
version of rppr::

  taxit reroot --rppr ~/local/bin/rppr -p my_refpkg


rollback
--------

.. literalinclude:: _helptext/rollback.txt


Examples:

Update the author on my_refpkg, then revert the change::

  taxit update --metadata 'author=Boris the mad baboon'
  taxit rollback my_refpkg

Roll back the last 3 operations on my_refpkg::

  taxit rollback -n 3 my_refpkg

rollforward
-----------

.. literalinclude:: _helptext/rollforward.txt


Examples:

Roll back the last operation on my_refpkg, then restore it::

  taxit rollback my_refpkg
  taxit rollforward my_refpkg

Roll forward the last 3 rollbacks on my_refpkg::

  taxit rollforward -n 3 my_refpkg

rp (resolve path)
-----------------

.. literalinclude:: _helptext/rp.txt

strip
-----

.. literalinclude:: _helptext/strip.txt

Examples:

Perform an update::

  taxit update my_refpkg hilda=file1

After this, file1 is still in the refpkg, but not referred to except
by the rollback information::

  taxit update my_refpkg hilda=file2

Now ``strip`` deletes file1, and the rollback and rollforward information::

  taxit strip my_refpkg


taxids
------

.. literalinclude:: _helptext/taxids.txt

Examples:

Look up two species and print their tax_ids to stdout, one per line::

  taxit taxids ncbi_database.db -n "Lactobacillus crispatus,Lactobacillus helveticus"

Read the species from some_names.txt and write their tax_ids to some_taxids.txt::

  taxit taxids ncbi_database.db -f some_names.txt -o some_taxids.txt

taxtable
--------

.. literalinclude:: _helptext/taxtable.txt

Examples:

Extract tax_ids 47770 and 33945 and all nodes connecting them to the root.::

  taxit taxtable taxonomy.db -t 47770,33945

The same as above, but write the output to subtax.csv instead of stdout::

  taxit taxtable taxonomy.db -t 47770,33945 -o subtax.csv

Extract the same tax_ids, plus the taxa specifies in taxnames.txt::

  taxit taxtable taxonomy.db -t 47770,33945 -n taxnames.txt -o taxonomy_from_both.csv

update
------

.. literalinclude:: _helptext/update.txt


update_taxids
-------------

.. literalinclude:: _helptext/update_taxids.txt


