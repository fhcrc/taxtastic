
exception Found_multiply

open MapsSets

(* remove exactly one element from l; fail if found multiply *)
let list_remove1 x l = 
  match List.partition ((=) x) l with
  | ([_],l') -> l'
  | ([],_) -> raise Not_found
  | (_,_) -> raise Found_multiply

let safe_add h k v = 
  assert(not (Hashtbl.mem h k));
  Hashtbl.add h k v


(* *** PTREE BUILDING *** *)

type edge = 
  (*        id    bl      connections *)
  | Pend of int * float * int list
  | Inte of int * float * int list * int list

  (* oh yes! *)
type ptree = (int, edge) Hashtbl.t

let of_stree bl_getter st = 
  let pt = Hashtbl.create (1+(Stree.max_id st)) in
  let add_edge id left righto = 
    let bl = bl_getter id in
    Hashtbl.add pt id
      (match righto with
      | None -> Pend(id, bl, left)
      | Some right -> Inte(id, bl, left, right))
  in
  let rec aux above_ids = function
    | Stree.Node(id, tL) ->
        let ids_of_tl = List.map Stree.top_id in
        add_edge id above_ids (Some (ids_of_tl tL));
        List.iter
          (fun (below, rest) -> aux (ids_of_tl rest) below) 
          (Base.pull_each_out tL)
    | Stree.Leaf id -> add_edge id above_ids None
  in
  aux [] st;
  Hashtbl.remove pt (Stree.top_id st); (* remove fake root leaf *)
  pt

let of_gtree gt = 
  of_stree 
    (fun i -> (IntMap.find i gt.Gtree.bark_map)#get_bl)
    gt.Gtree.stree

let of_string s = of_gtree (Newick.of_string s)


(* *** PTREE CHANGING *** *)

let autoedge id bl l1 l2 = 
  match (l1, l2) with
  | ([], []) -> assert_false
  | ([], l) | (l, []) -> Pend(id, bl, l)
  | _ -> Inte(id, bl, l1, l2)

(* not so pretty *)
let remove_id e rem_id = 
  match e with
  | Pend(id, bl, el) -> Pend(id, bl, list_remove1 el i)
  | Inte(id, bl, el1, el2) ->
      if List.mem rem_id el1 then autoedge id bl (list_remove1 rem_id el1) el2
      else autoedge id bl el1 (list_remove1 rem_id el2)

let sorted_list_eq l1 l2 = List.sort compare l1 = List.sort compare l2

(* assert that one_side is a (potentially resorted version of) l or r, and
 * return the one that it is not *)
let other_side one_side l r =
  if sorted_list_eq one_side l then r
  else if sorted_list_eq one_side r then l
  else assert(false)

(* replace the binding of k in h with its image under f *)
let hashtbl_map1 f h k = Hashtbl.replace h k (f (Hashtbl.find h k))



(* sound of rubber hitting road *)
let delete_pend pt del_id = 
  match Hashtbl.find pt del_id with
  | Inte(_,_,_,_) -> failwith "can't delete internal edge"
  | Pend(id, del_bl, el) ->
      assert(del_id == id);
      match el with
      | [] | [_] -> assert false 
      | [eid1; eid2] -> begin
        (* degree two-- heal the wound. *)
        match (Hashtbl.find pt eid1, Hashtbl.find pt eid2) with
        | (Pend(_,_,_),Pend(_,_,_)) -> assert(false)
        | (Inte(_,bl1,l1,r1), Inte(_,bl2,l2,r2)) -> 
        (* join two actual internal edges together. 
         * we arbitrarily pick eid1 for the id of the renewed edge. *)
            Hashtbl.replace pt eid1
              (Inte(eid1, bl1+.bl2, 
                other_side [del_id; eid2] l1 r1,
                other_side [del_id; eid1] l2 r2));
            Hashtbl.remove pt eid2
        | (Pend(_,bl1,l1), Inte(_,bl2,l2,r)) ->
            let new_node = 
              list_remove1 
other_side [del_id; eid2] l1 r
        (* we are deleting one side of a cherry. in this case, we extend the
         * branch length on the other pendant edge. *)
            Hashtbl.replace pt eid1
              (Pend(eid1, bl1+.del_id, );
        (* then we attach the new extended pendant directly to the edge leading
         * to the cherry *)
            Hashtbl.replace pt eid2
              (autoedge eid2 bl2 
            let other = other_side [del_id; eid1] l2 r in
        | (Inte(_,bl1,l1,r), Pend(_,bl2,l2)) ->
            Hashtbl.replace pt eid2
              (Pend(eid2, bl1+.bl2, other_side [del_id; eid1] l2 r);
            Hashtbl.remove pt eid1








(* we only want Pends in the set! *)
module OrderedPend = struct
  type t = edge
  let compare a b = match (a,b) with
  | (Pend(_, bla, _), Pend(_, blb, _)) -> compare bla blb
  | _ -> assert(false)
end

module PendSet = Set.Make(OrderedPend)

let pendset_of_pt pt = 
  Hashtbl.fold 
    (fun _ e s -> 
      match e with 
      | Pend(_,_,_) as pe -> PendSet.add pe s 
      | Inte(_,_,_,_) -> s)
    pt 
    PendSet.empty

let pl_of_hash h = Hashtbl.fold (fun k v l -> (k,v)::l) h []
