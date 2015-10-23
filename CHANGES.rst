==========================
 change log for taxtastic
==========================

0.5.4-dev
=========
 * new ``taxit taxid_classified`` that decides if a tax_id is primary and valid (True/False)
 * ``taxit update_taxids`` will halt on unknown tax_ids unless ``--unknowns FILE`` is specified
 * ``taxit update_taxids`` only requires a csv file with 'tax_id' column
 * ``taxit update_taxids`` takes an optional ``--name-column`` to assist in assigning tax_ids

0.5.4
=====

 * Add ``taxit taxtable --full`` for outputing all ranks in header.
 * Update subcommand help text
 * Generate Sphinx docs using help text emitted by subcommands (GH-70)

0.5.3
=====

 * Suppress warning when updating refpkg ``tree_stats`` file via ``taxit update``.

0.5.2
=====

 * Fix GH-63: "empirical_frequencies" now set to false when parsing FastTree AA statistics files
 * Close GH-64: "empirical_frequencies" is now available as a flag for PhyML statistics files
 * Fixed bug that prevented temporary files from being deleted

0.5.1
=====

 * Fix GH-62: "empirical_frequencies" was not set when parsing PhyML AA statistics files.

0.5.0
=====

 * Add ``.drop()`` ``.collapse()`` methods to ``taxtastic.taxtable.TaxNode``
 * Change ``is_classified`` column in taxonomy database: now does not mark
   below species as unclassified if the species-level classification is valid. [GH-59]
 * Add ``taxit composition`` - shows the taxonomic composition of a reference package at a given rank
 * Fix broken ``taxit lonelynodes``
 * Add ``taxit merge`` - Identifies tax_ids which have been merged, suggests new tax_ids.
 * Add ``taxit add_to_taxtable`` - adds nodes to a taxonomy [GH-60]
 * Fix support for newer versions of PhyML [GH-61]
 * Updates for compatibility with RAxML 7.7.2


0.4
===

 * 'names' table in the taxonomy database has a new column
   'is_classified' indicating whether 'tax_name' should be considered
   "classified".
 * Bugfix in ``taxit findcompany``
 * Support stdin as a source for ``taxit findcompany``
 * Taxonomy objects use NCBI ranks by default
 * Reference packages are created optionally (fixes creation of empty reference
   packages for commands like ``taxit info nonexisting.refpkg``)
 * Support zipped reference packages
 * Add a taxtable API: ``taxtastic.taxtable``
 * Remove some tests requiring a full taxonomy database
 * Rerooting reference packages on creation [GH-57]
 * More intelligent file name generation on clash
 * Deprecate the default ``create=True`` in ``taxtastic.refpkg.Refpkg``
 * Some PEP8 fixes


0.3.2
=====

 * version number contains abbreviated git sha identifying the commit.
 * Initial release to PyPI
 * Add findcompany subcommand
 * Add refpkg_intersection subcommand
 * Remove some obsolete components
 * Check required fields in seqinfo file [GH-46]
 * Add option to build taxtable from seqinfo file [GH-55]
 * Add subcommand to update taxids [GH-56]
 * Support FastTree AA and DNA log files
 * Fix rank order bug (infraorder was below parvorder)
 * Documentation updates
