# Docker images for taxtastic

See https://github.com/fhcrc/taxtastic

Images must be identified using a specific tag, for example:

```
% docker pull nghoffman/taxtastic:0.7.0
% docker run -it --rm nghoffman/taxtastic:0.7.0 taxit --version
taxit v0.7.0
% docker run -it --rm -v $(pwd):/working nghoffman/taxtastic:0.7.0 taxit info /working/urogenital-2017-03-30_named-1.0.refpkg
number of sequences: 8049
package components
aln_fasta
aln_sto
phylo_model
profile
seq_info
taxonomy
tree
tree_stats
```

