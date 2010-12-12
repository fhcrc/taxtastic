
open Pd

let pl_of_hash h = Hashtbl.fold (fun k v l -> (k,v)::l) h []

let x = list_remove1 3 [4;3;5]

let gt = Newick.of_string "(x:15,(a:3,b:4):3):1"
let pt = of_gtree gt
let x = pl_of_hash pt
let ps = pendset_of_pt pt
let x = PendSet.elements ps

let st = Gtree.get_stree (Newick.of_string "((x,y),(a,b),(c,d))")
let pt = of_stree (fun _ -> 1.) st
let x = pl_of_hash pt

(*
 * [(5, Inte (1., [6; 2], [3; 4])); (4, Pend (4, 1., [5; 3]));
 *  (3, Pend (3, 1., [5; 4])); (2, Inte (1., [6; 5], [0; 1]));
 *  (1, Pend (1, 1., [2; 0])); (0, Pend (0, 1., [2; 1]))] 
 *  *)

let ps = pendset_of_pt pt
let x = PendSet.elements ps

let pt = of_file "COG0001.auto1.fast.tre"
(* let out = perform pt 1e-2 *)
