open Core

module Term = struct
  module X = struct
    type t = {
      id : string ;
      name : string ;
      namespace : string ;
      is_a : string list ;
    }
    [@@deriving show, sexp]

    let compare t1 t2 = String.compare t1.id t2.id
  end

  include X

  let of_obo_term (ot:Obo.term) = {id = ot.id; name = ot.name ; namespace = ot.namespace ; is_a = ot.is_a}

  let get_parents t = t.is_a

  let get_name t = t.name

  module Set = Set.Make(X)

end


type t = {
  terms : Term.t array ;
  id_index : int String.Map.t ;
  is_a : int list array ;
}
type ontology = t

type domain = Biological_process | Molecular_function | Cellular_component
[@@deriving sexp]

let find_term o id =
  match Map.find o.id_index id with
  | None -> None
  | Some i -> Some o.terms.(i)

let of_obo (obo:Obo.t) ns =
  let exception Unknown_parent of string * string in
  let nso =
    match ns with
    | Biological_process -> "biological_process"
    | Cellular_component -> "cellular_component"
    | Molecular_function -> "molecular_function"
  in
  let terms =
    List.filter obo ~f:(fun ot -> String.equal ot.namespace nso) (*take only this namespace*)
    |> List.map ~f:(fun ot -> Term.of_obo_term ot) (* transform Obo.term in Term.t*)
    |> List.dedup_and_sort ~compare:Term.compare (* remove duplicate terms if any *)
  in
  let id_index =
    List.mapi terms ~f:(fun i (tt : Term.t) -> (tt.id, i))
    |> String.Map.of_alist_exn (*create string.map, throw exception but safe because of earlier dedup and Term comparison function*)
  in
  try
    let terms = Array.of_list terms in
    let is_a = Array.map terms ~f:(fun (t : Term.t) ->
        List.map t.is_a ~f:(fun parent_id ->
            match Map.find id_index parent_id with
            | None -> raise (Unknown_parent (t.id, parent_id))
            | Some i -> i
          )
      )
    in
    Ok { terms ; id_index ; is_a }
  with Unknown_parent (term, parent_term) ->
    Error (sprintf "Term %s has a unknown parent named %s" term parent_term)

let expand_term_list o tl =
  let rec add_term_to_closure ts t =
    if Set.mem ts t then ts
    else (
        let new_set = Set.add ts t in
        List.fold t.Term.is_a ~init:new_set ~f:(fun term_set t_id ->
            match  find_term o t_id with
            | None -> assert false (*should not happen if obo file correct*)
            | Some tt -> add_term_to_closure term_set tt
          )
      )
  in
  List.fold tl ~init:Term.Set.empty ~f:(fun ts t -> add_term_to_closure ts t)
  |> Set.to_list


let expand_id_list o il =
  let tl = List.filter_map il ~f:(fun i -> find_term o i) in
  let et = expand_term_list o tl in
  List.map et ~f:(fun t -> t.id)


let filter_terms o sl =
  List.filter sl ~f:(fun x -> not (Option.is_none (find_term o x)))

let term_names o =
  Array.map o.terms ~f:(fun t -> t.id, Term.get_name t)
  |> Array.to_list
  |> String.Map.of_alist_exn

module PKey = struct
  type t = int
  type pkey = t
  let compare = Int.compare

  let get_term o i = o.terms.(i)

  let find_term_by_id o id =
    Map.find o.id_index id

  let is_a_transitive_closure o l =
    let rec traverse acc x =
      if Set.mem acc x then acc
      else
        List.fold o.is_a.(x) ~init:(Set.add acc x) ~f:traverse
    in
    let leaves = Int.Set.of_list l in
    Set.fold leaves ~init:Int.Set.empty ~f:traverse
    |> Set.to_list

  module Map = Int.Map

  module Table = struct
    type 'a t = {
      values : 'a array ;
      terms : Term.t array ;
    }
    let create (ontology : ontology) v = {
      values = Array.create ~len:(Array.length ontology.terms) v ;
      terms = ontology.terms ;
    }
    let get t i = t.values.(i)
    let set t i v = t.values.(i) <- v
    let fold t ~init ~f =
      Array.fold2_exn t.terms t.values ~init ~f:(fun acc term value -> f ~term ~value acc)
  end
end
