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

  module Set = Set.Make(X)

end


type t = {
  id_to_term : Term.t String.Map.t ; 
}

let find_term o id = String.Map.find o.id_to_term id 

let has_absent_parents sm =
  let absent_parents_for_term sm t =
    List.exists t.Term.is_a ~f:(fun x -> Option.is_none (String.Map.find sm x))
  in
  String.Map.exists sm ~f:(fun t -> absent_parents_for_term sm t)

let of_obo (obo:Obo.t) ns =
  let nso =
    match ns with
    | `Biological_process -> "biological_process"
    | `Cellular_component -> "cellular_component"
    | `Molecular_function -> "molecular_function"
  in
  let itt =
    List.filter obo ~f:(fun ot -> String.equal ot.namespace nso) (*take only this namespace*)
    |> List.map ~f:(fun ot -> Term.of_obo_term ot) (* transform Obo.term in Term.t*)
    |> List.dedup_and_sort ~compare:Term.compare (* remove duplicate terms if any *)
    |> List.map ~f:(fun (tt:Term.t) -> (tt.id, tt)) (*transform term list in list of tuples (id, term)*)
    |> String.Map.of_alist_exn (*create string.map, throw exception but safe because of earlier dedup and Term comparison function*)
  in
  if has_absent_parents itt then Error "obo file is incomplete!" else  Ok {id_to_term = itt} 


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
