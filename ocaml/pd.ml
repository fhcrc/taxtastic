
exception Found_multiply
exception Other_side of int list * int list * int list 
exception Side_without of int * int list * int list 
exception Int_not_found of int
exception Not_implemented of string

open MapsSets

(* remove exactly one element from l; fail if found multiply *)
let list_remove1 x l = 
  match List.partition ((=) x) l with
  | ([_],l') -> l'
  | ([],_) -> raise (Int_not_found x)
  | (_,_) -> raise Found_multiply

let safe_add h k v = 
  assert(not (Hashtbl.mem h k));
  Hashtbl.add h k v


(* *** PTREE BUILDING *** *)

type edge = 
  (*        id    bl      connections *)
  | Pend of int * float * int list
  | Inte of       float * int list * int list

  (* oh yes! *)
type ptree = (int, edge) Hashtbl.t

let sorted_list_eq l1 l2 = List.sort compare l1 = List.sort compare l2

(* assert that one_side is a (potentially resorted version of) l or r, and
 * return the one that it is not *)
let get_other_side one_side l r =
  if sorted_list_eq one_side l then r
  else if sorted_list_eq one_side r then l
  else raise (Other_side (one_side, l, r))

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

let of_stree bl_getter st = 
  let pt = Hashtbl.create (1+(Stree.max_id st)) in
  let add_edge id left righto = 
    let bl = bl_getter id in
    Hashtbl.add pt id
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
  let root_build (to_build, rest) = aux (List.map Stree.top_id rest) to_build in
  let () = 
    match st with
    | Stree.Leaf _ -> ()
    | Stree.Node(_, [t1; t2]) -> begin
      List.iter root_build (Base.pull_each_out [t1; t2]);
      let (eid1, eid2) = (Stree.top_id t1, Stree.top_id t2) in
      match ((eid1,Hashtbl.find pt eid1), (eid2,Hashtbl.find pt eid2)) with
      | ((id1, Inte(bl1,l1,r1)), (id2, Inte(bl2,l2,r2))) -> 
          (* remove degree two node at root *)
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


(* uuuuuugly! *)
exception Found_inte of int
let find_internal pt = 
  try
    Hashtbl.iter 
      (fun i -> function | Inte(_,_,_) -> raise (Found_inte i) | _ -> ())
      pt;
    assert false
  with
  | Found_inte i -> i

let to_gtree pt = 
  (* start with an index bigger than anything in pt *)
  let count = ref (Hashtbl.fold (fun i _ -> max i) pt 0) 
  and m = ref IntMap.empty
  in
  let add_bl i bl = 
    m := 
      IntMap.add i 
        (new Newick_bark.newick_bark (`Of_bl_name_boot(Some bl, None, None)))
        !m
  in
  let rec aux ~bad our_id = 
    match Hashtbl.find pt our_id with
    | Inte(bl,l,r) ->
        incr count;
        add_bl (!count) bl;
        let our_side = get_side_without bad l r in
        Stree.Node(!count, List.map (aux ~bad:our_id) our_side)
    | Pend(id,bl,_) -> add_bl id bl; Stree.Leaf id
  in
  let start_edge = find_internal pt in
  match Hashtbl.find pt start_edge with
  | Inte(bl,l,r) -> 
      let stl = aux ~bad:(List.hd r) start_edge
      and str = aux ~bad:(List.hd l) start_edge
      in
      add_bl (Stree.top_id stl) (bl/.2.);
      add_bl (Stree.top_id str) (bl/.2.);
      incr count;
      Gtree.gtree (Stree.Node(!count, [stl;str])) !m
  | Pend(_,_,_) -> assert(false)

(* *** PTREE CHANGING *** *)


(* replace the binding of k in h with its image under f *)
let hashtbl_map1 f h k = Hashtbl.replace h k (f (Hashtbl.find h k))

let freplace f h id = Hashtbl.replace h id (f (Hashtbl.find h id))

(* <sound of rubber hitting road> *)
let delete_pend pt del_id = 
  match Hashtbl.find pt del_id with
  | Inte(_,_,_) -> failwith "can't delete internal edge"
  | Pend(_, _, el) -> 
      Hashtbl.remove pt del_id;
      (match el with
      | [] | [_] -> assert false 
      | [eid1; eid2] -> begin
        (* degree two-- heal the wound. *)
        match ((eid1,Hashtbl.find pt eid1), (eid2,Hashtbl.find pt eid2)) with
        | ((_,Pend(_,_,_)),(_,Pend(_,_,_))) -> 
            raise (Not_implemented "can't make trees smaller than three leaves")
        | ((id1, Inte(bl1,l1,r1)), (id2, Inte(bl2,l2,r2))) -> 
        (* join two actual internal edges together. *)
            Hashtbl.replace pt id1
              (Inte(bl1+.bl2, 
                get_other_side [del_id; id2] l1 r1,
                get_other_side [del_id; id1] l2 r2));
            Hashtbl.remove pt id2
        | ((pid, Pend(orig_id,bl1,l1)), (iid, Inte(bl2,l2,r)))
        | ((iid, Inte(bl2,l2,r)), (pid, Pend(orig_id,bl1,l1))) ->
        (* we are deleting one side of a cherry. in this case, we extend the
         * branch length on the other pendant edge. *)
            assert(sorted_list_eq l1 [iid; del_id]);
            Hashtbl.replace pt iid
              (Pend(orig_id, bl1+.bl2, get_other_side [del_id; pid] l2 r));
            Hashtbl.remove pt pid;
      end
      | eidl -> 
        (* degree greater than two: 
        * delete del_id from the edge lists of the other nodes *)
          let check_rem1 l = 
            let out = List.filter ((<>) del_id) l in
        (* make sure that internal nodes still have degree greater than two *)
            assert(1 < List.length out);
            out
          in
          List.iter 
            (freplace
              (function
                | Pend(id,bl,l) -> Pend(id,bl, check_rem1 l)
                | Inte(bl,l,r) -> Inte(bl, check_rem1 l, check_rem1 r))
              pt)
            eidl
      )

      
(* we only want Pends in the set! *)
module OrderedPend = struct
  type t = int * edge
  (* order first by bl, then by compare *)
  let compare a b = match (a,b) with
  | ((_,Pend(_, bla, _)), (_,Pend(_, blb, _))) -> begin 
    match compare bla blb with 0 -> compare a b | x -> x 
    end
  | _ -> assert(false)
end

module PendSet = Set.Make(OrderedPend)

let pendset_of_pt pt = 
  Hashtbl.fold 
    (fun id e s -> 
      match e with 
      | Pend(_,_,_) as pe -> PendSet.add (id,pe) s 
      | Inte(_,_,_) -> s)
    pt 
    PendSet.empty

let pl_of_hash h = Hashtbl.fold (fun k v l -> (k,v)::l) h []

let perform orig_pt stopping_bl = 
  let pt = Hashtbl.copy orig_pt in
  let rec aux accu s = 
    match PendSet.min_elt s with
    | (id, Pend(orig_id,bl,_) as p) ->
        if bl > stopping_bl then accu
        else begin
          Printf.printf "%d %g\n" orig_id bl;
          try
            delete_pend pt id; 
            aux ((orig_id,bl,pl_of_hash pt)::accu) (PendSet.remove p s)
          with
          | Not_found -> print_endline "lolo"; (orig_id,bl,pl_of_hash pt)::accu
        end
    | (_, Inte(_,_,_)) -> assert false
  in
  List.rev (aux [] (pendset_of_pt pt))


