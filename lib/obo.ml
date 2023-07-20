open Core

type term = {
  id : string ;
  name : string ;
  namespace : string ;
  is_a : string list ;
}
[@@deriving show]

type t = term list
[@@deriving show]

type stanza_type =
  | Term
  | Typedef
  | Instance

type parsed_line =
  | Pair of (string * string)
  | Header of stanza_type
  | Empty_line

let remove_comment s =
  Stdlib.String.trim (
    match String.lsplit2 s ~on:'!' with
    | Some (x, _) -> x
    | None -> s
  )

let parsed_line_of_string s =
  match remove_comment s with
  | "" -> Ok Empty_line
  | "[Term]" -> Ok (Header Term)
  | "[Typedef]" -> Ok (Header Typedef)
  | "[Instance]" -> Ok (Header Instance)
  | uc  ->
    match String.lsplit2 uc ~on:':' with
    | Some (x, y) -> Ok (Pair (x, Stdlib.String.trim y))
    | None -> Error "not an obo file"

let parse_lines l =
  Result.all (List.map l ~f:parsed_line_of_string)

let get_field l k =
  List.filter_map l ~f:(fun (x, y) -> if String.equal x k then Some y else None)

(*
let get_at_most_one_field l k =
  match get_field l k with
  | [] -> Ok None
  | [x] -> Ok (Some x)
  | _ -> Error "more than one"
*)

let get_exactly_one_field l k =
  match get_field l k with
  | [x] -> Ok x
  | [] -> Error "nothing found"
  | _ -> Error "more than one"

let term_of_list l =
  let open Let_syntax.Result in
  let* id = get_exactly_one_field l "id" in
  let* name = get_exactly_one_field l "name" in
  let+ namespace = get_exactly_one_field l "namespace" in
  let  is_a = get_field l "is_a" in
  {id ; name ; namespace ; is_a }

let break_between_lines l1 l2 =
  match (l1, l2) with
  | (_, Header _) -> true
  | _ -> false

let process_pair p =
  match p with
  | Pair (x, y) -> Some (x,y)
  | _ -> None

let process_block_of_lines l =
  match l with
  | Header Term :: t -> Some (term_of_list (List.filter_map t ~f:process_pair))
  | _  -> None

let process_list l =
  List.group l ~break:break_between_lines
  |> List.filter_map ~f:process_block_of_lines
  |> Result.all

let of_obo_file path =
  let open Let_syntax.Result in
  let* raw_lines = Utils.read_lines path in
  let* parsed_lines = parse_lines raw_lines in
  process_list parsed_lines

let filter_namespace o ~namespace =
  List.filter o ~f:(fun t -> String.equal t.namespace namespace)
