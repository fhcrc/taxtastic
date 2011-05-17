---------
TAXTASTIC
---------

We love it, but what is it?
---------------------------

Taxtastic software written in python used to build and maintain `reference package`_ s-- i.e. collections of reference trees, reference alignments, profiles, and associated taxonomic information.

A script named ``taxit`` provides a command line interface::

  % ./taxit -h
  usage: taxit.pyc [-h] [-V] [-v] [-q]
		   {taxtable,convexify,create,check,badgraph,reroot,help} ...

  Creation, validation, and modification of reference packages for use with
  `pplacer` and related software.

  positional arguments:
    {taxtable,convexify,create,check,badgraph,reroot,help}
      help                Detailed help for actions using `help <action>`
      check               Not yet implemented
      create              Creates a reference package
      taxtable            Creates a CSV file describing lineages for a set of
			  taxa
      convexify           Finds discordance in a reference package
      reroot              Reroots a reference package
      badgraph            Generates a csv file representing the ratio of leaves
			  to edges in a tree.

  optional arguments:
    -h, --help            show this help message and exit
    -V, --version         Print the version number and exit
    -v, --verbose         Increase verbosity of screen output (eg, -v is
			  verbose, -vv more so)
    -q, --quiet           Suppress output

 
Here's a quick synopsis of how it's used. 
Say we would like to assemble a reference package for the gene *rpoB*.
We have a reference alignment (here rpoB.fasta) and an ML tree built from the reference alignment (here RAxML_result.rpoB.PROTGAMMAWAGF).
We also have a whitespace-delimited list of taxon ids (here rpoB.ids_only) for the sequences in our reference alignment.
We also have a CSV file mapping sequence names to taxon ids (here rpoB.seq_info.csv).
Specifically, it contains two columns, labeled "seqname" and "tax_id" giving that mapping.

We can make a taxonomically-informed reference set in just two steps. 
First we generate the taxonomy::

  taxit taxtable -t rpoB.ids_only -o rpoB.taxonomy

When performed for the first time, this downloads the NCBI taxonomy and builds a SQLite3 database from it. 
That takes a while, but it only needs to be performed once (after that you can point to the database with -d).
Then::

  taxit create \
        --package-name rpoB.refpkg \
        --locus rpoB \
	--author "your name <you@some.address.edu>" \
	--package-version 1.0 \
	--tree-stats RAxML_info.rpoB.PROTGAMMAWAGF \
        --tree-file RAxML_result.rpsB.PROTGAMMAWAGF \
        --aln-fasta rpoB.aln \
        --taxonomy rpoB.taxonomy \
        --seq-info rpoB.seq_info.csv 


Now you can run pplacer to get taxonomically annotated placements, and visualizations combining both taxonomic and phylogenetic placements with two more steps::

  pplacer -c rpoB.refpkg sample.fasta
  guppy fat -c rpoB.refpkg sample.place


.. Targets ..
.. _reference package: http://github.com/fhcrc/taxtastic/wiki/refpkg
