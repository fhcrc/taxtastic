===========
 TAXTASTIC
===========

Taxtastic is a python package used to build and maintain reference
packages-- i.e. collections of reference trees, reference alignments,
profiles, and associated taxonomic information.

Installing
==========

``taxtastic`` requires Python 2.7.  The simplest method of installing
is using `pip <http://pip-installer.org>`_::

    pip2 install taxtastic

We strongly recommend installation into a virtualenv.


sqlite3
-------

Taxtastic uses recursive common table expressions to query the
taxonomy database, which requires that the Python ``sqlite3`` module
is built against sqlite3 library version of 3.8.3 or higher
(http://www.sqlite.org/releaselog/3_8_3.html). You can check the
version like this::

  python -c 'import sqlite3; print sqlite3.sqlite_version'

``python setup.py`` will exit with an error if the sqlite3 library
dependency is not met. On older systems, it is possible to replace the
builtin ``sqlite3`` module by installing ``pysqlite2`` with updated
sqlite3 libraries using a provided script (assuming an active
virtualenv)::

  dev/install_pysqlite.sh

After the script completes, confirm that ``pysqlite2`` was installed::

  python -c 'from pysqlite2 import dbapi2; print dbapi2.sqlite_version'

At this point, taxtastic may be installed as described above.


We love it, but what is it?
===========================

* quickstart_
* `full documentation`_

A script named ``taxit`` provides a command line interface::

   % ./taxit  --help

   usage: taxit [-h] [-V] [-v] [-q]
                {help,add_nodes,add_to_taxtable,check,composition,count_taxids,create,findcompany,info,lonelynodes,merge,merge_taxtables,new_database,refpkg_intersection,reroot,rollback,rollforward,rp,strip,taxids,taxtable,update,update_taxids}
                ...

   Creation, validation, and modification of reference packages for use with
   `pplacer` and related software.

   positional arguments:
     {help,add_nodes,add_to_taxtable,check,composition,count_taxids,create,findcompany,info,lonelynodes,merge,merge_taxtables,new_database,refpkg_intersection,reroot,rollback,rollforward,rp,strip,taxids,taxtable,update,update_taxids}
       help                Detailed help for actions using `help <action>`
       add_nodes           Add new nodes to a database
       add_to_taxtable     Add nodes to an existing taxtable csv
       check               Validate a reference package
       composition         Show taxonomic composition of a reference package
       count_taxids        Count tax_id appearances in a taxtable lineage
       create              Create a reference package
       findcompany         Find company for lonely nodes
       info                Show information about reference packages.
       lonelynodes         Extracts tax ids of all lonely nodes in a taxtable
       merge               Identify merged taxids and provide replacements
       merge_taxtables     Combine two or more taxtables
       new_database        Download NCBI taxonomy and create a database
       refpkg_intersection
                           Find the intersection of a taxtable and a refpkg's
                           taxonomy.
       reroot              Taxonomically reroots a reference package
       rollback            Undo an operation performed on a refpkg
       rollforward         Restore a change to a refpkg immediately after being
                           reverted
       rp                  Resolve path; get the path to a file in the reference
                           package
       strip               Remove rollback and rollforward information from a
                           refpkg
       taxids              Convert a list of taxonomic names into a list of
                           tax_ids
       taxtable            Create a tabular representation of taxonomic lineages
       update              Add or modify files or metadata in a refpkg
       update_taxids       Update obsolete tax_ids

   optional arguments:
     -h, --help            show this help message and exit
     -V, --version         Print the version number and exit
     -v, --verbose         Increase verbosity of screen output (eg, -v is
                           verbose, -vv more so)
     -q, --quiet           Suppress output

.. Targets ..
.. _quickstart: http://fhcrc.github.com/taxtastic/quickstart.html
.. _full documentation: http://fhcrc.github.com/taxtastic/index.html
