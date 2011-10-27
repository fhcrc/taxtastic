I have my files, just tell me how to put them in a refpkg!
==========================================================

You're impatient, you're a phylogeneticist, and you just need to get a refpkg put together so you can start running ``pplacer``.  You already have your sequences aligned and a tree built from them.  Okay, I get it.  First, here's a checklist of the files you need (and the arbitrary filenames I'll be calling them by for the reset of this page).  You'll have to take care of any format conversions yourself:

* ``tax_ids.txt``: A list of all the tax_ids that the sequences in your refpkg use, one per line in a text file.  Blank lines and Python style comments in this file will be ignored, and repeated entries are only considered once.  
* ``cmalign_profile``: If you used ``cmalign`` to produce your multiple, the profile you used.
* ``tree.nwk``: The tree you built from your sequences, in Newick format.
* ``tree_stats.txt``: The log file from RAxML or FastTree when it ran.
* ``seqs.fasta`` and ``seqs.sto``: The multiply aligned sequences in both FASTA and Stockholm formats.
* ``seq_info.csv``: A CSV file describing each of the aligned sequences.  It begins with one line giving the field names::

      "seqname","accession","tax_id","species_name","is_type"
  
  This is followed by one line for each sequence in your multiple alignment.  ``seqname`` should match the ID of the sequence in the FASTA file.  ``accession`` is a database reference for the sequence, which can be the same as ``seqname``.  In our work, ``seqname`` is an RDP accession number and ``accession`` is the NCBI accession number corresponding to that RDP entry.  ``tax_id`` is the entry in the taxonomy this sequence is mapped to, and ``species_name`` is the name associated with that entry.  ``is_type`` indicates whether this sequence is from a typestrain or not (should be ``TRUE`` or ``FALSE``).  If you aren't dealing with typestrains, just set it to ``FALSE`` for all your sequences.

There, have all that assembled?  Now it's a few steps to construct a refpkg.

First we have to turn the list of tax_ids into a usable taxonomy.  You will need a copy of the NCBI taxonomy in a format taxtastic understands in order to extract the minimal taxonomy containing the tax_ids in ``tax_ids.txt``.  Run::

    taxit new_database -d taxonomy.db

This will download the NCBI taxonomy into the current directory and load it into ``taxonomy.db``.  This takes a while!  Plan to go get a cup of coffee.  Then run::

    taxit taxtable -d taxonomy.db -t tax_ids.txt -o taxa.csv

``taxa.csv`` is a CSV file containing the minimum subtaxonomy of ``taxonomy.db`` which contains all the entries in ``tax_ids.txt``.

Now we are ready to build the refpkg. Assuming you're building a refpkg of the locus ``locus_name`` (usually something like ``16s_rRNA``) and you want to put the resulting refpkg in ``my_locus.refpkg``, run the command::

    taxit create -l locus_name -P my_locus.refpkg \
        --taxonomy taxa.csv \
        --aln-fasta seqs.fasta \
        --aln-sto seqs.sto \
        --seq-info seq_info.csv \
        --profile cmalign_profile \
        --tree-stats tree_stats.txt \
        --tree-file tree.nwk

You may also want to add some metadata fields using the options ``--author``, ``--description``, and ``--package-version``, but these are optional.  The above command is the bare minimum you need for a refpkg to be valid for use with ``pplacer``.
