==========================
 change log for taxtastic
==========================

0.5.0-dev
=========

 * Add ``.drop()`` ``.collapse()`` methods to ``taxtastic.taxtable.TaxNode``
 * Change ``is_classified`` column in taxonomy database: now does not mark
   below species as unclassified if the species-level classification is valid. [GH-59]
 * Add ``taxit composition`` - shows the taxonomic composition of a reference package at a given rank
 * Fix broken ``taxit lonelynodes``
 * Add ``taxit merge`` - Identifies tax_ids which have been merged, suggests new tax_ids.


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
