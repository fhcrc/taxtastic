
# Taxonomic algorithms for sequence selection

## Introduction
- In order to do phylogenetic placement, we need a reference tree.
- In order to make such a tree, we need to choose a set of taxa.
- A good reference set is one which is correct and allows us to accurately classify query sequences
- We would like this tree to be
  - robust and contain only high quality, correctly labeled sequences
  - very diverse, to aid in high-resolution phylogenetic placement
  - have specific taxa of interest
  - of a reasonable size
- Potential objective functions
  - for "spread," we want to minimize the maximial branch length of a placed sequence
    - would be cool to do a study of all of the given methods, choosing between them based on this criterion
  - phylogenetic diversity
  - maximal internal branch length
  - percentage of phylogenetic diversity
  - for reliability [note that excluding difficult to place taxa may bias our taxon set because these taxa may have an "interesting" evol history]
  - maximal bootstrap/clade probability values
  - taxa have low delta values
  - tree is an agreement subtree of trees in the credible set

### Ref set selection goal
We don't want an enormous tree with everything in it because

* don't want disproportionate representation
* rampant mislabeling and poor quality sequences
* in the future this problem will get even worse
* why include something that you will never see?

### A priori selection criteria
* have as few placements as possible on the branches proximal to the given one
* shortest average pendant branch length for placements
* as many placements as possible end up on the leaves
* "widest" set of taxa, with deepest MRCA
* remove mislabeled taxa

### Potentially implementable selection criteria
* max PD
  * is maximum expected pairwise PD the same as max total PD?
* for "widest", we could try to maximize the expected depth of the MRCA of the selected taxa
  * problem is that this doesn't select the "deepest" member of a given clade
* we could maximize expected pairwise node distance
  * this would bias things towards heavily sampled taxa, as those taxa would have lots of nodes

### Usage scenarios
* what are all of the lineages associated with a given taxid
* merge locally generated sequences with previous reference set
* converting one set of names to another using a synonym table indicating the preferred name
* update taxonomic table to a new taxonomy (could involve reclassification of tax ids)
* remove mislabeled sequences (OK to call ocaml code? call RAxML to recompute tree?)
* given a big tree and alignment, cut it down to size to have at most X sequences per taxon
* add "off target" sequences as well as contaminant sequences


Research questions
* Say we take a given model of tree shape. Can we prove that a given means of selecting taxa minimizes the expected branch length to the "reference" tree?
* If we have a number of trees in the credible set, can we take an agreement subtree which has high PD, or some other desirable criterion?




Attic
=====
- split compatibility idea:
  - take all of the splits of the trees in the credible set, each of which are equipped with a distance
  - warm up idea: take the largest (in terms of weight sum) compatible subset of splits for our tree topology
  - can select a set of taxa with the same idea as follows:
    - selection of a taxon subset gives a projection on splits, and throw away the splits which project to something trivial
    - then select splits to maximize the total weight, such that the projection of the collection of splits is compatible

