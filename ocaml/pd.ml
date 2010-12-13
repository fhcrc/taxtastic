(* taxtastic v0.1. Copyright (C) 2009-2010  Frederick A Matsen. 
 * This file is part of taxtastic. Taxtastic is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version. Taxtastic is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details. You should have received a copy of the GNU General Public License along with taxtastic.  If not, see <http://www.gnu.org/licenses/>.  *)

open MapsSets
open Ptree

let freplace f h id = Hashtbl.replace h id (f (Hashtbl.find h id))

type idbl = {id : int; bl : float}
      
(* we only want Pends in the set! *)
module OrderedIdbl = struct
  type t = idbl
  (* order first by bl, then by id *)
  let compare a b = 
    match compare a.bl b.bl with 0 -> compare a.id b.id | x -> x 
end

module IdblSet = Set.Make(OrderedIdbl)

let idblset_of_pt pt = 
  Hashtbl.fold 
    (fun id e s -> 
      match e with 
      | Pend(_,bl,_) -> IdblSet.add {id=id;bl=bl} s 
      | Inte(_,_,_) -> s)
    pt 
    IdblSet.empty

(* <sound of rubber hitting road> *)
let delete_pend pt del_id idbls = 
  match Hashtbl.find pt del_id with
  | Inte(_,_,_) -> failwith "can't delete internal edge"
  | Pend(_, _, eidl) -> 
      Hashtbl.remove pt del_id;
      (match eidl with
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
            Hashtbl.remove pt id2;
            idbls
        | ((pid, Pend(orig_id,bl1,l1)), (iid, Inte(bl2,l2,r)))
        | ((iid, Inte(bl2,l2,r)), (pid, Pend(orig_id,bl1,l1))) ->
        (* we are deleting one side of a cherry. in this case, we extend the
         * branch length on the other pendant edge. *)
            assert(sorted_list_eq l1 [iid; del_id]);
            Hashtbl.replace pt iid
              (Pend(orig_id, bl1+.bl2, get_other_side [del_id; pid] l2 r));
            Hashtbl.remove pt pid;
            IdblSet.add 
              {id=iid; bl=bl1+.bl2} 
              (IdblSet.remove {id=pid; bl=bl1} idbls)
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
            eidl;
          idbls
      )

let pl_of_hash h = Hashtbl.fold (fun k v l -> (k,v)::l) h []

let perform orig_pt stopping_bl = 
  let pt = Hashtbl.copy orig_pt in
  let rec aux accu s = 
    let m = IdblSet.min_elt s in
    if m.bl > stopping_bl then accu
    else match Hashtbl.find pt m.id with
    | Pend(orig_id, bl, _) ->
        assert(bl = m.bl);
        let new_s = delete_pend pt m.id (IdblSet.remove m s) in 
        aux ((orig_id,m.bl,to_stree pt)::accu) new_s
    | Inte(_,_,_) -> assert false
  in
  List.rev (aux [] (idblset_of_pt pt))


