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

type domain = Biological_process | Molecular_function | Cellular_component

let find_term o id =
  match String.Map.find o.id_index id with
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
            match String.Map.find id_index parent_id with
            | None -> raise (Unknown_parent (t.id, parent_id))
            | Some i -> i
          )
      )
    in
    Ok { terms ; id_index ; is_a }
  with Unknown_parent (term, parent_term) ->
    Error (sprintf "Term %s has a unknown parent named %s" term parent_term)

let define_domain domain =
  match domain with
  | "biological_process" -> Ok Biological_process
  | "cellular_component" -> Ok Cellular_component
  | "molecular_function" -> Ok Molecular_function
  | _ -> Error (Printf.sprintf "Unknown GO domain %s, exiting.\n" domain)

let expand_term_list o tl =
  let rec add_term_to_closure ts t =
    if Term.Set.mem ts t then ts
    else (
        let new_set = Term.Set.add ts t in
        List.fold t.is_a ~init:new_set ~f:(fun term_set t_id ->
            match  find_term o t_id with
            | None -> assert false (*should not happen if obo file correct*)
            | Some tt -> add_term_to_closure term_set tt
          )
      )
  in
  List.fold tl ~init:Term.Set.empty ~f:(fun ts t -> add_term_to_closure ts t)
  |> Term.Set.to_list


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
