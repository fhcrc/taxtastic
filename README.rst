---------
TAXTASTIC
---------

Installing
----------

``taxtastic`` requires Python 2.7.  The simplest method of installing is using `pip <http://pip-installer.org>`_::

    pip2 install taxtastic

If you don't have pip, try::

    easy_install taxtastic

Either of these commands will install taxtastic along with required dependencies.

We love it, but what is it?
---------------------------

Taxtastic is software written in python used to build and maintain reference packages-- i.e. collections of reference trees, reference alignments, profiles, and associated taxonomic information.

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
