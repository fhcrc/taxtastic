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
  usage: taxit.pyc [-h] [-V] [-v] [-q]

                   {info,rollback,help,create,update,new_database,reroot,taxids,strip,add_nodes,rollforward,check,taxtable}
                   ...

  Creation, validation, and modification of reference packages for use with
  `pplacer` and related software.

  positional arguments:
    {info,rollback,help,create,update,new_database,reroot,taxids,strip,add_nodes,rollforward,check,taxtable}
      help                Detailed help for actions using `help <action>`
      add_nodes           Add new nodes to a database containing a taxonomy.
      check               taxtastic/subcommands/check.py Run a series of deeper
                          checks on a RefPkg. This subcommand is a wrapper
                          around the Refpkg method is_ill_formed.
      info                Show information about reference packages.
      create              Creates a reference package
      new_database        Creates a CSV file describing lineages for a set of
                          taxa
      reroot              Taxonomically reroots a reference package
      update              Adds or updates files or metdata in a refpkg. The
                          update subcommand takes a refpkg to operate on, then a
                          series of changes to make, expressed as key=file. So
                          to add a file ../otherdir/boris under the key 'meep'
                          and abcd under the key 'hilda' in a refpkg 'my-
                          refpkg', you would run $ taxit update my-refpkg
                          meep=../otherdir/boris hilda=abcd If a file already
                          exists under a given key, it is overwritten. Passing
                          taxit update the --metadata option makes it update the
                          metadata instead of files. For example, to set the
                          author field to "Genghis Khan" and the version to
                          0.4.3, run $ taxit update --metadata "author=Genghis
                          Khan" version=0.4.3
      taxids              Look up a set of tax_ids from taxonomic names
      taxtable            Creates a CSV file describing lineages for a set of
                          taxa
      strip               Removes all rollback and rollforward information and
                          files not attached to the current state from a refpkg.
                          $ taxit strip my-refpkg
      rollback            Rollback a refpkg to undo the previous command. $
                          taxit rollback my-refpkg You can also specify -n # to
                          specify the number of operations to roll back
                          (defaults to 1), as in $ taxit rollback -n 3 my-refpkg
      rollforward         Rollforward a rolled back command on a refpkg. $ taxit
                          rollforward my-refpkg You can also specify -n # to
                          specify the number of operations to roll forward
                          (defaults to 1), as in $ taxit rollforward -n 3 my-
                          refpkg

  optional arguments:
    -h, --help            show this help message and exit
    -V, --version         Print the version number and exit
    -v, --verbose         Increase verbosity of screen output (eg, -v is
                          verbose, -vv more so)
    -q, --quiet           Suppress output

  usage: taxit.pyc [-h] [-V] [-v] [-q]


.. Targets ..
.. _quickstart: http://fhcrc.github.com/taxtastic/quickstart.html
.. _full documentation: http://fhcrc.github.com/taxtastic/index.html
