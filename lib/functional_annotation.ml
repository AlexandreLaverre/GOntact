open Core
    
type t = {
  gene_symbol_to_go : string list String.Map.t ;
  gene_id_to_go : string list String.Map.t ;
  go_to_gene_id : string list String.Map.t ;
  go_to_gene_symbol : string list String.Map.t ;  
}

let extract_terms ga is =
  let d, k = 
    match is with
    | `Id id ->  ga.gene_id_to_go, id
    | `Symbol sym -> ga.gene_symbol_to_go, sym
  in
  String.Map.find d k 

let extract_terms_exn  ga is =
  let d, k = 
    match is with
    | `Id id ->  ga.gene_id_to_go, id
    | `Symbol sym -> ga.gene_symbol_to_go, sym
  in
  String.Map.find_exn d k 

let extract_genes ga ~go_id:go i =
  let d = match i with
    | `Id -> ga.go_to_gene_id
    | `Symbol -> ga.go_to_gene_symbol
  in
  String.Map.find d go

let of_gaf_and_ontology (gaf:Gaf.t) (ont:Ontology.t) =
  let only_this_namespace = List.filter gaf ~f:(fun ga -> not (Option.is_none (Ontology.find_term ont ga.go_id))) in
  let gaf_without_no = List.filter only_this_namespace ~f:(fun ga -> not (String.is_prefix ga.qualifier ~prefix:"NOT")) in
  let no_empty_symbols = List.filter gaf_without_no ~f:(fun ga -> not (String.equal ga.gene_symbol "")) in
  let make_dict fg f compare =
    let m = String.Map.of_alist_multi (List.map fg ~f) in
    String.Map.map m ~f:(List.dedup_and_sort ~compare) 
  in
  let gene_symbol_to_go = make_dict no_empty_symbols (fun ga -> ga.gene_symbol, ga.go_id) String.compare in
  let gene_id_to_go = make_dict gaf_without_no (fun ga -> (ga.gene_id, ga.go_id)) String.compare in
  let go_to_gene_id = make_dict gaf_without_no (fun ga -> (ga.go_id,  ga.gene_id)) String.compare in
  let go_to_gene_symbol = make_dict no_empty_symbols (fun ga -> (ga.go_id,  ga.gene_symbol)) String.compare in
  {gene_symbol_to_go ; gene_id_to_go ; go_to_gene_symbol ; go_to_gene_id }

let show fa d =
  match d with
    | `Gene_id_to_go -> String.Map.to_alist fa.gene_id_to_go |> [%show: (string * string list) list] 
    | `Gene_symbol_to_go -> String.Map.to_alist fa.gene_symbol_to_go |>  [%show: (string * string list) list] 
    | `Go_to_gene_id -> String.Map.to_alist fa.go_to_gene_id |>  [%show: (string * string list) list]
    | `Go_to_gene_symbol -> String.Map.to_alist fa.go_to_gene_symbol |>  [%show: (string * string list) list]


let reverse_dict gtg =
  let al = String.Map.to_alist gtg in (*transform string.map in list of tuples (key, val)*)
  let flatten_list (k, l) =
    List.fold l ~init:[] ~f:(fun ll x -> (x, k)::ll)
  in
  let reval = List.fold al ~init:[] ~f:(fun ll x -> List.append (flatten_list x) ll) in 
  String.Map.of_alist_multi reval 

let propagate_annotations fa o =
  let filter_and_extend dict =
    let filtered = String.Map.map dict ~f:(fun l -> Ontology.filter_terms o l) in  (* for each gene, extract GO ids that belong to this ontology *) 
    String.Map.map filtered ~f:(fun l -> Ontology.expand_id_list o l) (*expand id list *)
  in
  let extended_gig = filter_and_extend fa.gene_id_to_go in
  let extended_gsg = filter_and_extend fa.gene_symbol_to_go in
  let extended_go_to_gene_id = reverse_dict extended_gig in (*create expanded go id to gene dictionaries *)
  let extended_go_to_gene_symbol = reverse_dict extended_gsg in
  {gene_symbol_to_go = extended_gsg; gene_id_to_go = extended_gig; go_to_gene_symbol = extended_go_to_gene_symbol; go_to_gene_id = extended_go_to_gene_id}


let write_annotations fa path =
  let s_to_go = fa.gene_symbol_to_go in 
  let output = Out_channel.create ~append:false path in
  let write_entry k l = 
     List.iter l ~f:(Printf.fprintf output "%s\t%s\n" k) 
  in
  Out_channel.output_string output "GeneSymbol\tGOID\n" ; 
  String.Map.iteri s_to_go ~f:(fun ~key:k ~data:l -> write_entry k l) ; 
  Out_channel.close output 
  
let gene_symbols fa =
  String.Map.keys fa.gene_symbol_to_go
  |> String.Set.of_list

let go_list_of_gene_symbol fa gs =
  let symbol2go = fa.gene_symbol_to_go in
  match String.Map.find symbol2go gs with
  | Some l -> l
  | None -> [] 
            
