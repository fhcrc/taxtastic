
exception Found_multiply

open MapsSets

(* remove exactly one element from l; fail if found multiply *)
let list_remove_single x l = 
  match List.partition ((=) x) l with
  | ([_],l') -> l'
  | ([],_) -> raise Not_found
  | (_,_) -> raise Found_multiply

let safe_add h k v = 
  assert(not (Hashtbl.mem h k));
  Hashtbl.add h k v

type edge = 
  (*        id    bl      connections *)
  | Pend of int * float * int list
  | Inte of int * float * int list * int list

  (* oh yes! *)
type ptree = (int, edge) Hashtbl.t

let ptree_of_stree bl_getter st = 
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

let ptree_of_gtree gt = 
  ptree_of_stree 
    (fun i -> (IntMap.find i gt.Gtree.bark_map)#get_bl)
    gt.Gtree.stree

let pl_of_hash h = Hashtbl.fold (fun k v l -> (k,v)::l) h []
