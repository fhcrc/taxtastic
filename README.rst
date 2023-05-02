===========
 TAXTASTIC
===========

Taxtastic is a python package used to build and maintain reference
packages, i.e. collections of reference trees, reference alignments,
profiles, and associated taxonomic information.

.. image:: https://travis-ci.org/fhcrc/taxtastic.svg?branch=master
    :target: https://travis-ci.org/fhcrc/taxtastic

We love it, but what is it?
===========================

* quickstart_
* `full documentation`_

A script named ``taxit`` provides a command line interface::

  % taxit  --help
  usage: taxit [-h] [-V] [-v] [-q]
	       {help,add_nodes,add_to_taxtable,check,composition,create,extract_nodes,findcompany,get_lineage,info,lonelynodes,new_database,refpkg_intersection,reroot,rollback,rollforward,rp,strip,taxids,taxtable,update,update_taxids}
	       ...

  Creation, validation, and modification of reference packages for use with
  `pplacer` and related software.

  positional arguments:
    {help,add_nodes,add_to_taxtable,check,composition,create,extract_nodes,findcompany,get_lineage,info,lonelynodes,new_database,refpkg_intersection,reroot,rollback,rollforward,rp,strip,taxids,taxtable,update,update_taxids}
      help                Detailed help for actions using `help <action>`
      add_nodes           Add nodes and names to a database
      add_to_taxtable     Add nodes to an existing taxtable csv
      check               Validate a reference package
      composition         Show taxonomic composition of a reference package
      create              Create a reference package
      extract_nodes       Extract nodes from a given source in yaml format
      findcompany         Find company for lonely nodes
      get_lineage         Calculate the taxonomic lineage of a taxid
      info                Show information about reference packages.
      lonelynodes         Extracts tax ids of all lonely nodes in a taxtable
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
      taxids              Convert a list of taxonomic names into a recursive
			  list of species
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
.. _quickstart: https://fhcrc.github.io/taxtastic/quickstart.html
.. _full documentation: https://fhcrc.github.io/taxtastic/index.html


Installation
============

``taxtastic`` requires Python 3.8+. The simplest method of installing
is using `pip <http://pip-installer.org>`_::

  pip install taxtastic

We strongly recommend installation into a virtualenv. Instructions for
installing the ``taxtastic`` package and the ``taxit`` command line
entry point in a virtualenv are as follows::

  python3 -m venv taxtastic-env
  source taxtastic-env/bin/activate
  pip install -U pip
  pip install taxtastic

If you prefer to install from the git repository::

  git clone https://github.com/fhcrc/taxtastic.git
  cd taxtastic
  python3 -m venv taxtastic-env
  source taxtastic-env/bin/activate
  pip install .

Finally, ``taxit`` can be run from a docker image hosted from Docker
Hub (https://hub.docker.com/r/nghoffman/taxtastic/), for example::

  docker run --rm -it -v $(pwd):$(pwd) -w $(pwd) nghoffman/taxtastic:latest taxit -v new_database

sqlite3
-------

Taxtastic uses recursive common table expressions to query the
taxonomy database, which requires that the Python ``sqlite3`` module
is built against sqlite3 library version of 3.8.3 or higher
(http://www.sqlite.org/releaselog/3_8_3.html). You can check the
version like this::

  python3 -c 'import sqlite3; print(sqlite3.sqlite_version)'

A note on databases
===================

This project supports both sqlite3 and postgresql as database
backends. For most applications, we recommend sqlite3: some operations
(particularly initial database creation) are much faster using sqlite3
due to the details of how postgresql enforces database constraints (we
may try to optimize this in the future - in theory, postgresql can be
made to be at least as fast). If you do want to use postgresql, note
that some of the queries consume a lot of memory, and the default
configuration tends to be memory constrained (and this *really* slows
things down). On a reasonably new mac laptop, we found that the
optimizations suggested here
(http://big-elephants.com/2012-12/tuning-postgres-on-macos/) do the
trick.
