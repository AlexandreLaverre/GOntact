open Core
(*
type go_association = {
  go_id : string ;
  qualifier : string ;
}
*)
    
type gene = {
  id : string ;
  symbol : string ;
}
[@@deriving ord]

(*
type gene_annotation = {
  gene : gene ;
  go_associations : go_association list ;
}
*)

type t = {
  gene_symbol_to_go : string list String.Map.t ;
  gene_id_to_go : string list String.Map.t ;
  go_to_gene : gene list String.Map.t ;  
}

let extract_terms ga is =
  let d, k = 
    match is with
    | `Id id ->  ga.gene_id_to_go, id
    | `Symbol sym -> ga.gene_symbol_to_go, sym
  in
  String.Map.find d k 

let extract_genes ga ~go_id:go i =
  let d = ga.go_to_gene in 
  let proj =
    match i with
    | `Id -> fun gene -> gene.id
    | `Symbol -> fun gene -> gene.symbol
  in 
  Option.map (String.Map.find d go) ~f:(List.map ~f:proj)


let from_gaf (gaf:Gaf.t) =
  let filtered_gaf = List.filter gaf ~f:(fun ga -> not (String.is_prefix ga.qualifier ~prefix:"NOT")) in 
  let make_dict f compare =
    let m = String.Map.of_alist_multi (List.map filtered_gaf ~f) in
    String.Map.map m ~f:(List.dedup_and_sort ~compare) 
  in
  let gene_symbol_to_go = make_dict (fun ga -> ga.gene_symbol, ga.go_id) String.compare in
  let gene_id_to_go = make_dict (fun ga -> (ga.gene_id, ga.go_id)) String.compare in
  let go_to_gene = make_dict (fun ga -> (ga.go_id, {id = ga.gene_id ; symbol = ga.gene_symbol})) compare_gene in
  {gene_symbol_to_go ; gene_id_to_go ; go_to_gene }
