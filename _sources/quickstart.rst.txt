Quickstart
==========


Minimal example
---------------

Say we have the following:

* ``seqs.fasta``: The multiply aligned reference sequences in FASTA format.
* ``tree.nwk``: The tree built from reference alignment, in Newick format.
* ``tree_stats.txt``: The log file from FastTree or the statistics file from RAxML/phyml.

Build a reference package ``my.refpkg`` for a locus ``locus_name`` (e.g. ``16s_rRNA``) in a single command as follows::

    taxit create -l 16s_rRNA -P my.refpkg \
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

  This is followed by one line for each sequence in your multiple alignment. Only ``seqname`` and ``tax_id`` are required. ``seqname`` should match the ID of the sequence in the FASTA file and must be unique. ``accession`` provides an additional identifier for the sequence, which can (but need not) be the same as ``seqname``.  In our work, ``seqname`` is often an RDP accession number and ``accession`` is the NCBI accession number corresponding to that RDP entry (note that an NCBI accession number may not uniquely identify a sequence, as in tha case of a genomic sequence with multiple 16S rRNA gene loci).``tax_id`` is the entry in the taxonomy this sequence is mapped to, and ``species_name`` is the name associated with that entry. ``is_type`` indicates whether this sequence is from a type strain or not (should be ``TRUE`` or ``FALSE``).  If you aren't dealing with typestrains, just set it to ``FALSE`` for all your sequences.

Now there are a few additional steps to construct a taxonomically-annotated refpkg.

First we have to represent the tax_ids as a taxonomy.  You will need a copy of the NCBI taxonomy in a format taxtastic understands in order to extract the minimal taxonomy containing the tax_ids in ``tax_ids.txt``.  Run::

    taxit new_database taxonomy.db

This will download zipped files containing the NCBI taxonomy into the current directory and load the data into an sqlite databse named ``taxonomy.db``.  This takes a while!  Plan to go get a cup of coffee.  Then run::

    taxit taxtable taxonomy.db -f tax_ids.txt -o taxa.csv

The output, ``taxa.csv``, is a CSV file containing the minimum subtaxonomy of ``taxonomy.db`` which contains all the entries in ``tax_ids.txt``.

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

where ``align_profile`` is a profile that can be used to align query sequences using HMMER3 or Infernal, and ``seqs.sto`` is the reference alignment in Stockholm format. This is handy if you plan to use ``cmalign`` (in the Infernal suite of tools) or ``hmmalign`` (in the HMMER3 suite of tools) to perform your alignments.
