
open Pd

let x = list_remove_single 3 [4;3;5]

let gt = Newick.of_string "(a:3.,b:4.):1"

let pt = of_gtree gt

let x = pl_of_hash pt
