
exception Found_multiply
exception Other_side of int list * int list * int list 
exception Int_not_found of int

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
  aux [] st;
  Hashtbl.remove pt (Stree.top_id st); (* remove fake root leaf *)
  pt

let of_gtree gt = 
  of_stree 
    (fun i -> (IntMap.find i gt.Gtree.bark_map)#get_bl)
    gt.Gtree.stree

let of_string s = of_gtree (Newick.of_string s)
let of_file s = of_gtree (Newick.of_file s)



(* *** PTREE CHANGING *** *)

let sorted_list_eq l1 l2 = List.sort compare l1 = List.sort compare l2

(* assert that one_side is a (potentially resorted version of) l or r, and
 * return the one that it is not *)
let other_side one_side l r =
  if sorted_list_eq one_side l then r
  else if sorted_list_eq one_side r then l
  else raise (Other_side (one_side, l, r))

(* replace the binding of k in h with its image under f *)
let hashtbl_map1 f h k = Hashtbl.replace h k (f (Hashtbl.find h k))

let freplace f h id = Hashtbl.replace h id (f (Hashtbl.find h id))

(* sound of rubber hitting road *)
let delete_pend pt del_id = 
  match Hashtbl.find pt del_id with
  | Inte(_,_,_) -> failwith "can't delete internal edge"
  | Pend(_, _, el) -> 
      Hashtbl.remove pt del_id;
      (match el with
      | [] -> assert false
      | [_] -> assert false 
      | [eid1; eid2] -> begin
        (* degree two-- heal the wound. *)
        match ((eid1,Hashtbl.find pt eid1), (eid2,Hashtbl.find pt eid2)) with
        | ((_,Pend(_,_,_)),(_,Pend(_,_,_))) -> assert false
        | ((id1, Inte(bl1,l1,r1)), (id2, Inte(bl2,l2,r2))) -> 
        (* join two actual internal edges together. *)
            Hashtbl.replace pt id1
              (Inte(bl1+.bl2, 
                other_side [del_id; id2] l1 r1,
                other_side [del_id; id1] l2 r2));
            Hashtbl.remove pt id2
        | ((pid, Pend(orig_id,bl1,l1)), (iid, Inte(bl2,l2,r)))
        | ((iid, Inte(bl2,l2,r)), (pid, Pend(orig_id,bl1,l1))) ->
        (* we are deleting one side of a cherry. in this case, we extend the
         * branch length on the other pendant edge. *)
            assert(sorted_list_eq l1 [iid; del_id]);
            Hashtbl.replace pt iid
              (Pend(orig_id, bl1+.bl2, other_side [del_id; pid] l2 r));
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
  let compare a b = match (a,b) with
  | ((_,Pend(_, bla, _)), (_,Pend(_, blb, _))) -> compare bla blb
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

let perform orig_pt stopping_bl = 
  let pt = Hashtbl.copy orig_pt in
  let rec aux accu s = 
    match PendSet.min_elt s with
    | (id, Pend(orig_id,bl,_) as p) ->
        if bl > stopping_bl then accu
        else begin
          Printf.printf "%d %g\n" orig_id bl;
          delete_pend pt id; 
          aux ((orig_id,bl)::accu) (PendSet.remove p s)
        end
    | (_, Inte(_,_,_)) -> assert false
  in
  List.rev (aux [] (pendset_of_pt pt))


