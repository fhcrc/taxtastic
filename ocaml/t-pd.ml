
open Pd

let pl_of_hash h = Hashtbl.fold (fun k v l -> (k,v)::l) h []

let x = list_remove1 3 [4;3;5]

let gt = Newick.of_string "(x:15,(a:3,b:4):3):1"
let pt = of_gtree gt
let x = pl_of_hash pt
let ps = pendset_of_pt pt
let x = PendSet.elements ps

let test s = 
  let st = Gtree.get_stree (Newick.of_string s) in
  let pt = of_stree (fun _ -> 1.) st in
  pl_of_hash pt

let x = test "((x,y),(a,b))"
let x = test "((x,y),(a,b),(c,d))"

let ps = pendset_of_pt pt
let x = PendSet.elements ps

let pt = of_file "COG0001.auto1.fast.tre"
(* let out = perform pt 1e-2 *)
