
open Pd

let pl_of_hash h = Hashtbl.fold (fun k v l -> (k,v)::l) h []

let x = list_remove1 3 [4;3;5]

let unit_bl_of_s s =
  let st = Gtree.get_stree (Newick.of_string s) in
  of_stree (fun _ -> 1.) st

let test s = 
  let pt = unit_bl_of_s s in
  let ps = pendset_of_pt pt in
  (pl_of_hash pt, PendSet.elements ps)

let x = test "((x,y),(a,b))"
let x = test "((x,y),(a,b),(c,d))"

let four = of_string "((x:1,y:2):4,(a:9,b:9):9,(c:9,d:9):9):9";;
let x = PendSet.elements (pendset_of_pt four);;
let x = pl_of_hash four;;
let x = perform four 7.;;

(* let pt = of_file "COG0001.auto1.fast.tre" *)
(* let out = perform pt 1e-2 *)
