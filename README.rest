---------
TAXTASTIC
---------

We love it, but what is it?
---------------------------

Taxtastic is a series of python scripts to build and maintain `reference package`_ s-- i.e. collections of reference trees, reference alignments, profiles, and associated taxonomic information.
We've only just begun.

Here's a quick synopsis of how it's used. 
Say we would like to assemble a reference package for the gene *rpoB*.
We have a reference alignment (here rpoB.fasta) and an ML tree built from the reference alignment (here RAxML_result.rpoB.PROTGAMMAWAGF).
We also have a whitespace-delimited list of taxon ids (here rpoB.ids_only) for the sequences in our reference alignment.
We also have a CSV file mapping sequence names to taxon ids (here rpoB.seq_info.csv).
Specifically, it contains two columns, labeled "seqname" and "tax_id" giving that mapping.

We can make a taxonomically-informed reference set in just two steps. 
First we generate the taxonomy::

  taxtable.py -t rpoB.ids_only > rpoB.taxonomy

When performed for the first time, this downloads the NCBI taxonomy and builds a SQLite3 database from it. 
That takes a while, but it only needs to be performed once (after that you can point to the database with -d).
Then::

  taxomatic.py create \
        -s RAxML_info.rpoB.PROTGAMMAWAGF \
        -t RAxML_result.rpsB.PROTGAMMAWAGF \
        -f rpoB.aln \
        -T rpoB.taxonomy \
        -i rpoB.seq_info.csv \
        -P rpoB.refpkg 

Now you can run pplacer to get taxonomically annotated placements, and visualizations combining both taxonomic and phylogenetic placements with two more steps::

  pplacer -c rpoB.refpkg sample.fasta
  guppy fat -c rpoB.refpkg sample.place



.. Targets ..
.. _reference package: http://github.com/fhcrc/taxtastic/wiki/refpkg
