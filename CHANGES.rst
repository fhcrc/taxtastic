==========================
 change log for taxtastic
==========================

0.4
===

 * 'names' table in the taxonomy database has a new column
   'is_classified' indicating whether 'tax_name' should be considered
   "classified".


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
