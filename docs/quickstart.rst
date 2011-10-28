I have my files, just tell me how to put them in a refpkg!
==========================================================

You're impatient, you're a phylogeneticist, and you just need to get a refpkg put together so you can start running ``pplacer``.  You already have your sequences aligned and a tree built from them.  Okay, I get it.  First, here's a checklist of the files you need (and the arbitrary filenames I'll be calling them by for the reset of this page).  You'll have to take care of any format conversions yourself:

* ``cmalign_profile``: If you used ``cmalign`` to produce your multiple, the profile you used.
* ``tree.nwk``: The tree you built from your sequences, in Newick format.
* ``tree_stats.txt``: The log file from RAxML or FastTree when it ran.
* ``seqs.fasta`` and ``seqs.sto``: The multiply aligned sequences in both FASTA and Stockholm formats.

There, have all that assembled?  Now we shove it in a refpkg.

Assuming you're building a refpkg of the locus ``locus_name`` (usually something like ``16s_rRNA``) and you want to put the resulting refpkg in ``my_locus.refpkg``, run the command::

    taxit create -l locus_name -P my_locus.refpkg \
        --aln-fasta seqs.fasta \
        --aln-sto seqs.sto \
        --profile cmalign_profile \
        --tree-stats tree_stats.txt \
        --tree-file tree.nwk

You may also want to add some metadata fields using the options ``--author``, ``--description``, and ``--package-version``, but these are optional.  The above command is the bare minimum you need for a refpkg to be valid for use with ``pplacer``.
