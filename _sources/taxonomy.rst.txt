Manipulating the NCBI taxonomy
==============================

The first thing we must do is download a copy of the NCBI taxonomy.  The taxit subcommand ``new_database`` does this for us and loads the taxonomy into an SQLite3 database that taxit can use.::

    $ taxit new_database taxonomy.db

The primary use for the taxonomy is extracting subsets of it for use in references packages.  For instance, if we want the minimal subtaxonomy that contains both Lactobacillus crispatus (tax_id 47770) and Enterococcus avium (tax_id 33945), we would run::

    $ taxit taxtable taxonomy.db -t 47770,33945 -o minimal_taxonomy.csv

The file ``minimal_taxonomy.csv`` contains the entries for the two species we specified, along with all nodes of the taxonomy which connect them to the root of the taxonomy.  The file's contents are::

    "tax_id","parent_id","rank","tax_name","root","below_root","superkingdom","phylum","class","order","family","genus","species"
    "1","1","root","root","1","","","","","","","",""
    "131567","1","below_root","cellular organisms","1","131567","","","","","","",""
    "2","131567","superkingdom","Bacteria","1","131567","2","","","","","",""
    "1239","2","phylum","Firmicutes","1","131567","2","1239","","","","",""
    "91061","1239","class","Bacilli","1","131567","2","1239","91061","","","",""
    "186826","91061","order","Lactobacillales","1","131567","2","1239","91061","186826","","",""
    "81852","186826","family","Enterococcaceae","1","131567","2","1239","91061","186826","81852","",""
    "33958","186826","family","Lactobacillaceae","1","131567","2","1239","91061","186826","33958","",""
    "1350","81852","genus","Enterococcus","1","131567","2","1239","91061","186826","81852","1350",""
    "1578","33958","genus","Lactobacillus","1","131567","2","1239","91061","186826","33958","1578",""
    "33945","1350","species","Enterococcus avium","1","131567","2","1239","91061","186826","81852","1350","33945"
    "47770","1578","species","Lactobacillus crispatus","1","131567","2","1239","91061","186826","33958","1578","47770"

Each entry is listed with the information about its node (tax_id, parent, rank, etc.), and a series of fields giving the lineage of taxids connecting this node to the root of the taxonomy.

For large numbers of tax_ids, it is awkard or impossible to pass them on the command line to ``taxit taxtable``.  In these cases, you can also give the ``-t`` option a filename containing taxids.  The file understands Python style comments, and extracts tax_ids from each line seperated by spaces, commas, or semicolons, as in::

    # This is an example file of tax_ids
    47770 # one per line works just fine
    33945 87541, 358;    1383 # or separations with spaces, commas, or semicolons

    # and empty lines are just fine
    47770 # repeating taxids is just fine; they will only be included one

If this were placed in ``taxids.txt``, we would extract the minimal taxonomy containing these taxids with::

    $ taxit taxtable taxonomy.db -f taxids.txt -o taxonomy_from_taxids.csv

Sometimes you will have the names of taxonomic nodes instead of tax_ids or some names and some tax_ids.  ``taxit taxtable`` also accepts names, passed in a file via the ``-n`` option.  For instance, if we had the tax_ids from the above file besides 47770 and 33945 as names instead,::

    Aerococcus christensenii # Python style comments work here as well
    Agrobacterium tumefaciens, Atopobium rimae # but only , and ; separate fields, not space
    Atopobium rimae # Again, repeated names don't matter

If this were placed in ``taxnames.txt`` and we passed the original tax_ids of 47770 and 33945 as well, we would extract the same subtaxonomy as above with::

    $ taxit taxtable taxonomy.db -t 47770,33945 -n taxnames.txt -o taxonomy_from_both.csv

In these examples we have used only species level tax_ids.  This is not necessary.  If you pass ``taxit taxtable`` a family level tax_id, it will add the family and all nodes connecting it to the root of the taxonomy.  However, it will not add any species, genera, or the like below the level of that family, so passing higher rank tax_ids in this manner is not terribly useful.

Occasionally you will want to augment the NCBI taxonomy with entries of your own.  The ``taxit add_nodes`` command takes an Excel 97 or CSV file, with a header line, and adds its entries to the database.  The new nodes must include fields for ``tax_id``, ``parent_id``, ``rank``, and ``tax_name``.  For instance, if we want to add a new species Lactobacillus borisii, to which we have assigned the arbitrary tax_id AA1 (the NCBI taxonomy doesn't use letters in its tax_ids, so a safe way to avoid collisions is to choose a couple letters as a prefix for your own additions), we would put it in the following form:::

    "tax_id","parent_id","rank","tax_name"
    "AA1","1578","species","Lactobacillus borisii"

where 1578 is the tax_id of the genus Lactobacillus.  If we put this text in ``new_taxa.csv``, we would add it with the command::

    $ taxit add_nodes taxonomy.db -N new_taxa.csv

Then we could refer to either Lactobacillus borisii or the tax_id AA1 when extracting subtaxonomies from ``taxonomy.db``.

