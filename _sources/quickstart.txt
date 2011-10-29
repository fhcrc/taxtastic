Quickstart
==========


Minimal example
---------------

Say we have the following:

* ``seqs.fasta``: The multiply aligned reference sequences in FASTA format.
* ``tree.nwk``: The tree built from reference alignment, in Newick format.
* ``tree_stats.txt``: The log file from FastTree or the statistics file from RAxML/phyml.

Build a reference package ``my.refpkg`` for a locus ``locus_name`` (e.g. ``16s_rRNA``) in a single command as follows::

    taxit create -l locus_name -P my.refpkg \
        --aln-fasta seqs.fasta \
        --tree-stats tree_stats.txt \
        --tree-file tree.nwk


Taxonomically equipped reference package example
------------------------------------------------

The fun really begins when we incorporate taxonomic information.
Say in addition to the above files you also have:

* ``tax_ids.txt``: A list of all the tax_ids that the sequences in your refpkg use, one per line in a text file.  Blank lines and Python style comments in this file will be ignored, and repeated entries are only considered once.
* ``seq_info.csv``: A CSV file describing each of the aligned sequences.  It begins with one line giving the field names::

      "seqname","accession","tax_id","species_name","is_type"

  This is followed by one line for each sequence in your multiple alignment.  ``seqname`` should match the ID of the sequence in the FASTA file.  ``accession`` is a database reference for the sequence, which can be the same as ``seqname``.  In our work, ``seqname`` is an RDP accession number and ``accession`` is the NCBI accession number corresponding to that RDP entry.  ``tax_id`` is the entry in the taxonomy this sequence is mapped to, and ``species_name`` is the name associated with that entry.  ``is_type`` indicates whether this sequence is from a typestrain or not (should be ``TRUE`` or ``FALSE``).  If you aren't dealing with typestrains, just set it to ``FALSE`` for all your sequences.

Now it's a few steps to construct a taxonomically-annotated refpkg.

First we have to turn the list of tax_ids into a usable taxonomy.  You will need a copy of the NCBI taxonomy in a format taxtastic understands in order to extract the minimal taxonomy containing the tax_ids in ``tax_ids.txt``.  Run::

    taxit new_database -d taxonomy.db

This will download the NCBI taxonomy into the current directory and load it into ``taxonomy.db``.  This takes a while!  Plan to go get a cup of coffee.  Then run::

    taxit taxtable -d taxonomy.db -t tax_ids.txt -o taxa.csv

``taxa.csv`` is a CSV file containing the minimum subtaxonomy of ``taxonomy.db`` which contains all the entries in ``tax_ids.txt``.

Now we are ready to build the refpkg::

    taxit create -l locus_name -P my.refpkg \
        --taxonomy taxa.csv \
        --aln-fasta seqs.fasta \
        --seq-info seq_info.csv \
        --tree-stats tree_stats.txt \
        --tree-file tree.nwk

You may also want to add some metadata fields using the options ``--author``, ``--description``, and ``--package-version``, but these are optional. In addition, you may package in profile information and a Stockholm alignment::

        --aln-sto seqs.sto \
        --profile align_profile

where ``align_profile`` is a profile that can be used to align query sequences using HMMER3 or Infernal, and ``seqs.sto`` is the reference alignment in Stockholm format.
