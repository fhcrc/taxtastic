(* taxtastic v0.1. Copyright (C) 2009-2010  Frederick A Matsen. 
 * This file is part of taxtastic. Taxtastic is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version. Taxtastic is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details. You should have received a copy of the GNU General Public License along with taxtastic.  If not, see <http://www.gnu.org/licenses/>.  *)

exception Other_side of int list * int list * int list 
exception Side_without of int * int list * int list 
exception Not_implemented of string

open MapsSets

type edge = 
  (*        id    bl      connections *)
  | Pend of int * float * int list
  | Inte of       float * int list * int list

(* oh yes! a mutable data structure, just for fun. *)
type ptree = (int, edge) Hashtbl.t


(* *** GENERALITIES *** *)

let safe_add h k v = 
  assert(not (Hashtbl.mem h k));
  Hashtbl.add h k v

let sorted_list_eq l1 l2 = List.sort compare l1 = List.sort compare l2

(* assert that one_side is a (potentially resorted version of) l or r, and
 * return the one that it is not *)
let get_other_side one_side l r =
  if sorted_list_eq one_side l then r
  else if sorted_list_eq one_side r then l
  else raise (Other_side (one_side, l, r))

(* same, but for a single element *)
let get_side_without x l r = 
  if List.mem x l then r
  else if List.mem x r then l
  else raise (Side_without (x, l, r))

let list_replace ~src ~dst = List.map (fun x -> if x = src then dst else x)

let edgel_replace ~src ~dst = function
  | Pend(id,bl,l) as p -> 
      if List.mem src l then Pend(id,bl, list_replace src dst l) else p
  | Inte(bl,l,r) as i -> 
      if not ((List.mem src l) || (List.mem src r)) then i
      else Inte(bl, list_replace src dst l, list_replace src dst r)


(* *** OF *** *)

(* bl_getter is a function which returns a bl given an id *)
let of_stree bl_getter st = 
  let pt = Hashtbl.create (1+(Stree.max_id st)) in
  (* righto is None if we have a pendant edge *)
  let add_edge id left righto = 
    let bl = bl_getter id in
    safe_add pt id
      (match righto with
      | None -> Pend(id, bl, left)
      | Some right -> Inte(bl, left, right))
  in
  let rec aux above_ids = function
    | Stree.Node(id, tL) ->
        let ids_of_tl = List.map Stree.top_id in
        add_edge id above_ids (Some (ids_of_tl tL));
        List.iter
          (fun (below, rest) -> aux (id::(ids_of_tl rest)) below) 
          (Base.pull_each_out tL)
    | Stree.Leaf id -> add_edge id above_ids None
  in
  (* the basic tree builder, which doesn't take tricky rooting into account *)
  let root_build (to_build, rest) = 
    aux (List.map Stree.top_id rest) to_build 
  in
  let () = 
    match st with
    | Stree.Leaf _ -> ()
    | Stree.Node(_, [t1; t2]) -> begin
      (* tree with degree two rootings require some special care *)
      (* first build the basic tree *)
      List.iter root_build (Base.pull_each_out [t1; t2]);
      let (id1, id2) = (Stree.top_id t1, Stree.top_id t2) in
      match (Hashtbl.find pt id1, Hashtbl.find pt id2) with
      | (Inte(bl1,l1,r1), Inte(bl2,l2,r2)) -> 
          let join1 = get_other_side [id2] l1 r1
          and join2 = get_other_side [id1] l2 r2
          in
          (* make new internal edge *)
          Hashtbl.replace pt id1 (Inte(bl1+.bl2, join1, join2));
          (* clean out old edge *)
          Hashtbl.remove pt id2;
          (* reconnect things to new edge *)
          List.iter 
            (fun id -> 
              Hashtbl.replace pt id 
                (edgel_replace ~src:id2 ~dst:id1 (Hashtbl.find pt id)))
            (join1 @ join2);
      | _ -> raise (Not_implemented "rooted on pendant edge")
      end
    | Stree.Node(_, tL) -> List.iter root_build (Base.pull_each_out tL);
  in
  pt

let of_gtree gt = 
  let get_bl i = 
    try (IntMap.find i gt.Gtree.bark_map)#get_bl with
    | Not_found -> failwith "tree is missing branch lengths"
  in
  of_stree get_bl gt.Gtree.stree

let of_string s = of_gtree (Newick.of_string s)
let of_file s = of_gtree (Newick.of_file s)


(* *** TO *** *)

(* uuuuuugly! *)
exception Found_inte of int
let find_internal pt = 
  try
    Hashtbl.iter 
      (fun i -> function 
        | Inte(_,_,_) -> raise (Found_inte i) 
        | _ -> ())
      pt;
    None
  with
  | Found_inte i -> Some i

let to_gtree pt = 
  (* start with maximum of indices of the ids *)
  let count = 
    ref (Hashtbl.fold 
      (fun _ e j -> match e with Pend(i,_,_) -> max i j | _ -> j) 
      pt 0) 
  and m = ref IntMap.empty
  in
  let add_bark i bl nameo = 
    m := 
      IntMap.add i 
        (new Newick_bark.newick_bark (`Of_bl_name_boot(Some bl, nameo, None)))
        !m
  in
  (* above is a neighboring edge index in the "up" direction in the resulting
   * rooting *)
  let rec aux ~above our_id = 
    match Hashtbl.find pt our_id with
    | Inte(bl,l,r) ->
        incr count;
        let node_id = !count in (* have to nail down count due to recursion *)
        add_bark node_id bl None;
        let our_side = get_side_without above l r in
        Stree.Node(node_id, List.map (aux ~above:our_id) our_side)
    | Pend(id,bl,_) -> add_bark id bl (Some (string_of_int id)); Stree.Leaf id
  in
  match find_internal pt with
  | Some start_edge -> begin
      match Hashtbl.find pt start_edge with
      | Inte(bl,l,r) -> 
          let stl = aux ~above:(List.hd r) start_edge
          and str = aux ~above:(List.hd l) start_edge
          in
          add_bark (Stree.top_id stl) (bl/.2.) None;
          add_bark (Stree.top_id str) (bl/.2.) None;
          incr count;
          Gtree.gtree (Stree.Node(!count, [stl;str])) !m
      | Pend(_,_,_) -> assert(false)
  end
  | None ->
      let tL = (* fix our mutables *)
        Hashtbl.fold
          (fun _ e l ->
            match e with 
            | Inte(_,_,_) -> assert false
            | Pend(id,bl,_) -> 
                add_bark id bl (Some (string_of_int id));
                (Stree.Leaf id)::l)
          pt
          []
      in
      Gtree.gtree (Stree.Node (1 + !count, tL)) !m

let to_stree pt = Gtree.get_stree (to_gtree pt)
