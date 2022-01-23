open Core

type gene_annotation = {
  gene_id : string ;
  gene_symbol : string ;  
  qualifier : string ;
  go_id :string ;
}
[@@deriving show]

type t = gene_annotation list
[@@deriving show]

(* column numbers for this GAF format
   let index_db_object_id = 0
   let index_db_object_symbol = 1
   let index_qualifier = 2
   let index_go_id = 3
 *)

let parse_line l =
  if String.is_prefix l ~prefix:"!" then Ok None
  else
    match String.split l ~on:'\t' with
    | [_ ;  gene_id ; gene_symbol ; qualifier ; go_id ; _ ; _ ; _ ; _ ; _ ; _ ; _ ; _ ; _ ; _ ; _ ; _ ] -> Ok (Some {gene_id; gene_symbol; qualifier; go_id})
    | _ -> Error "not a GAF file"
        
let from_file path = 
  In_channel.read_lines path
  |> List.map ~f:parse_line
  |> Result.all
  |> Result.map ~f:List.filter_opt
  
