---------
TAXTASTIC
---------

We love it, but what is it?
---------------------------

Taxtastic is software written in python used to build and maintain reference packages-- i.e. collections of reference trees, reference alignments, profiles, and associated taxonomic information.

* quickstart_
* `full documentation`_

A script named ``taxit`` provides a command line interface::

   % ./taxit -h


   usage: taxit [-h] [-V]

                {info,rollback,help,create,strip,taxids,new_database,check,reroot,refpkg_intersection,lonelynodes,update_taxids,rp,add_nodes,rollforward,update,findcompany,taxtable}
                ...

   Creation, validation, and modification of reference packages for use with
   `pplacer` and related software.

   positional arguments:
     {info,rollback,help,create,strip,taxids,new_database,check,reroot,refpkg_intersection,lonelynodes,update_taxids,rp,add_nodes,rollforward,update,findcompany,taxtable}
       help                Detailed help for actions using `help <action>`
       add_nodes           Add new nodes to a database containing a taxonomy.
       check               Run a series of deeper checks on a RefPkg.
       create              Creates a reference package
       findcompany         Find company for lonely nodes.
       info                Show information about reference packages.
       lonelynodes         Extracts tax ids of all lonely nodes in a taxtable.
       new_database        Create a database containing an entire taxonomy
       refpkg_intersection
                           Find the intersection of a taxtable and a refpkg's
                           taxonomy.
       reroot              Taxonomically reroots a reference package
       rollback            Rollback a refpkg to undo the previous command.
       rollforward         Rollforward a rolled back command on a refpkg.
       rp                  Resolve path; get the path to a file in the reference
                           package.
       strip               Removes all rollback and rollforward information and
                           files not attached to the current state from a refpkg.
       taxids              Look up a set of tax_ids from taxonomic names
       taxtable            Creates a CSV file describing lineages for a set of
                           taxa
       update              Adds or updates files or metdata in a refpkg.
       update_taxids       Update obsolete tax_ids in preparation for use in
                           taxtable. Takes sequence info

   optional arguments:
     -h, --help            show this help message and exit
     -V, --version         Print the version number and exit


.. Targets ..
.. _quickstart: http://fhcrc.github.com/taxtastic/quickstart.html
.. _full documentation: http://fhcrc.github.com/taxtastic/index.html
