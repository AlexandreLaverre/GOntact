open Core

(* invariant: list of sorted lists *)
type t = Union of Ontology.PKey.t list list

let of_sorted_lists_unsafe xs = Union xs

let sorted_list_union xs ys ~compare =
  let rec loop acc xs ys =
    match xs, ys with
    | [], ys -> List.rev_append acc ys
    | xs, [] -> List.rev_append acc xs
    | hx :: tx, hy :: ty ->
      match compare hx hy with
      |  0 -> loop (hx :: acc) tx ty
      | -1 -> loop (hx :: acc) tx ys
      |  1 -> loop (hy :: acc) xs ty
      | _ -> assert false
  in
  loop [] xs ys

let to_sorted_list (Union u) =
  match u with
  | [] -> []
  | [xs] -> xs
  | xs ->
    List.reduce_exn xs ~f:(sorted_list_union ~compare:Ontology.PKey.compare)

let union (Union u) (Union v) = Union (u @ v)

(* let rec sort_lists_by_head = function *)
(*   | [] -> [] *)
(*   | [] :: t -> sort_lists_by_head t *)
(*   | (h1 :: _ as l1) :: t -> *)
(*     match sort_lists_by_head t with *)
(*     | [] -> [ l1 ] *)
(*     | [] :: _ -> assert false *)
(*     | (h2 :: _ as l2) :: t as sorted_t -> *)
(*       if Ontology.PKey.compare h1 h2 <= 0 then l1 :: sorted_t *)
(*       else l2 :: (l1 :: t) *)

let rec traverse u ~f =
  match find_next u with
  | None -> ()
  | Some (e, u') ->
    f e ; traverse u' ~f
and find_next = function
  | [] -> None
  | [] :: rest -> find_next rest
  | (h :: t as l0) :: rest ->
    match find_better_candidate h rest with
    | None, rest' -> Some (h, t :: rest')
    | Some h', rest' -> Some (h', l0 :: rest')
and find_better_candidate x = function
  | [] -> None, []
  | [] :: rest -> find_better_candidate x rest
  | (h :: t as l) :: rest ->
    match Ontology.PKey.compare x h with
    | -1 -> (
        let ans, rest' = find_better_candidate x rest in
        ans, l :: rest'
      )
    | 0 -> (
        let ans, rest' = find_better_candidate x rest in
        match ans with
        | None -> None, t :: rest'
        | Some _ -> ans, l :: rest'
      )
    | 1 -> (
        let ans, rest' = find_better_candidate h rest in
        match ans with
        | None -> Some h, t :: rest'
        | Some _ -> ans, l :: rest'
      )
    | _ -> assert false

let rec traverse2 xs ys ~f =
  match xs, ys with
  | [], _ -> List.iter ys ~f
  | _, [] -> List.iter xs ~f
  | h_x :: t_x, h_y :: t_y ->
    match Ontology.PKey.compare h_x h_y with
    | -1 -> f h_x ; traverse2 t_x ys ~f
    |  0 -> f h_x ; traverse2 t_x t_y ~f
    |  1 -> f h_y ; traverse2 xs t_y ~f
    | _ -> assert false

let iter (Union u) ~f =
  match u with
  | [] -> ()
  | [xs] -> List.iter xs ~f
  | [xs ; ys] -> traverse2 xs ys ~f
  | _ -> traverse u ~f
