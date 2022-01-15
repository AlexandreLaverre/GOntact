open Core
    
type term = {
  id : string ;
  name : string ;
  namespace : string ;
  is_a : string list ;
}

type t = term list

type stanza_type =
  | Term
  | Typedef
  | Instance

type parsed_line =
  | Pair of (string * string)
  | Header of stanza_type
  | Empty_line

                   
let rec skip_until_term get_next_line =
  match get_next_line () with
  | None ->  `End_of_file
  | Some "[Term]" -> `Next_term
  | Some _ ->  skip_until_term get_next_line


let remove_comment s =
  match String.lsplit2 s ~on:'!' with
  | Some (x, _) -> x (* we should trim string *)
  | None -> s
    
let parsed_line_of_string s =
  match remove_comment s with
  | "" -> Ok Empty_line
    | "[Term]" -> Ok (Header Term)
    | "[Typedef]" -> Ok (Header Typedef)
    | "[Instance]" -> Ok (Header Instance)
    | _  ->
      match String.lsplit2 s ~on:':' with
      | Some (x, y) -> Ok (Pair (x, y))
      | None -> Error "not an obo file"
        
let rec parse_lines l =
  Result.all (List.map l ~f:parsed_line_of_string)

let get_field l k =
  List.filter_map l ~f:(fun (x, y) -> if String.equal x k then Some y else None)

let get_at_most_one_field l k =
  match get_field l k with
  | [] -> Ok None
  | [x] -> Ok (Some x)
  | _ -> Error "more than one"

let get_exactly_one_field l k =
  match get_field l k with
  | [x] -> Ok x
  | [] -> Error "nothing found"
  | _ -> Error "more than one"

  
let term_of_list l =
  let* id = get_exactly_one_field l "id"


    with  
      {id = ;  name = ; namespace = ; is_a = } 

(* let parse_term get_next_line =  *)
  
  
let parse get_next_line =
  match skip_until_term get_next_line with
  | `End_of_file -> []   (* file does not contain any GO term, empty list *)
  | `Next_term -> parse_term get_next_line
  | _ -> print_endline "weird"
    
    
(* let from_file path = *)
             
  
