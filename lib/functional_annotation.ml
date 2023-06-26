open Core

type t = {
  gene_symbol_to_go : Ontology.PKey.t list String.Map.t ;
  gene_id_to_go : Ontology.PKey.t list String.Map.t ;
  go_to_gene_id : string list Ontology.PKey.Map.t ;
  go_to_gene_symbol : string list Ontology.PKey.Map.t ;
  ontology : Ontology.t ;
}

let term_names_of_pkeys fa xs =
  List.map xs ~f:(fun pkey -> (Ontology.PKey.get_term fa.ontology pkey).id)

let extract_terms ga is =
  let d, k =
    match is with
    | `Id id ->  ga.gene_id_to_go, id
    | `Symbol sym -> ga.gene_symbol_to_go, sym
  in
  String.Map.find d k

let extract_terms_exn ga is =
  match extract_terms ga is with
  | None -> raise Stdlib.Not_found
  | Some xs -> xs

let extract_genes ga ~go_id:go i =
  let d = match i with
    | `Id -> ga.go_to_gene_id
    | `Symbol -> ga.go_to_gene_symbol
  in
  Ontology.PKey.find_term_by_id ga.ontology go
  |> Option.bind ~f:(Ontology.PKey.Map.find d)

let of_gaf_and_ontology (gaf:Gaf.t) (ont:Ontology.t) =
  let only_this_namespace = List.filter_map gaf ~f:(fun ga ->
      match Ontology.PKey.find_term_by_id ont ga.go_id with
      | None -> None
      | Some key -> Some (ga, key)
    ) in
  let gaf_without_no = List.filter only_this_namespace ~f:(fun (ga, _) -> not (String.is_prefix ga.qualifier ~prefix:"NOT")) in
  let no_empty_symbols = List.filter gaf_without_no ~f:(fun (ga, _) -> not (String.is_empty ga.gene_symbol)) in
  let make_pkey_dict xs proj =
    List.fold xs ~init:Ontology.PKey.Map.empty ~f:(fun acc x ->
        let key, data = proj x in
        Ontology.PKey.Map.add_multi acc ~key ~data
      )
    |> Ontology.PKey.Map.map ~f:(List.dedup_and_sort ~compare:String.compare)
  in
  let make_sm_dict xs proj =
    List.fold xs ~init:String.Map.empty ~f:(fun acc x ->
        let key, data = proj x in
        String.Map.add_multi acc ~key ~data
      )
    |> String.Map.map ~f:(List.dedup_and_sort ~compare:Ontology.PKey.compare)
  in
  let gene_symbol_to_go = make_sm_dict no_empty_symbols (fun (ga, k) -> ga.gene_symbol, k) in
  let gene_id_to_go = make_sm_dict gaf_without_no (fun (ga, k) -> (ga.gene_id, k)) in
  let go_to_gene_id = make_pkey_dict gaf_without_no (fun (ga, k) -> (k,  ga.gene_id)) in
  let go_to_gene_symbol = make_pkey_dict no_empty_symbols (fun (ga, k) -> (k,  ga.gene_symbol)) in
  { gene_symbol_to_go ; gene_id_to_go ; go_to_gene_symbol ; go_to_gene_id ; ontology = ont }

let show fa d =
  let get_id pkey = (Ontology.PKey.get_term fa.ontology pkey).id in
  match d with
  | `Gene_id_to_go ->
    String.Map.map fa.gene_id_to_go ~f:(term_names_of_pkeys fa)
    |> String.Map.to_alist
    |> [%show: (string * string list) list]
  | `Gene_symbol_to_go ->
    String.Map.map fa.gene_symbol_to_go ~f:(term_names_of_pkeys fa)
    |> String.Map.to_alist
    |>  [%show: (string * string list) list]
  | `Go_to_gene_id ->
    fa.go_to_gene_id
    |> Ontology.PKey.Map.to_alist
    |> List.map ~f:(fun (pkey, xs) -> get_id pkey, xs)
    |> [%show: (string * string list) list]
  | `Go_to_gene_symbol ->
    fa.go_to_gene_symbol
    |> Ontology.PKey.Map.to_alist
    |> List.map ~f:(fun (pkey, xs) -> get_id pkey, xs)
    |> [%show: (string * string list) list]


let propagate_annotations fa o =
  let reverse_dict gtg =
    String.Map.fold gtg ~init:Ontology.PKey.Map.empty ~f:(fun ~key ~data acc ->
        List.fold data ~init:acc ~f:(fun acc pkey ->
            Ontology.PKey.Map.add_multi acc ~key:pkey ~data:key
          )
      )
    |> Ontology.PKey.Map.map ~f:(List.dedup_and_sort ~compare:String.compare)
  in
  let filter_and_extend dict =
    String.Map.map dict ~f:(Ontology.PKey.is_a_transitive_closure o) (*expand id list *)
  in
  let extended_gig = filter_and_extend fa.gene_id_to_go in
  let extended_gsg = filter_and_extend fa.gene_symbol_to_go in
  let extended_go_to_gene_id = reverse_dict extended_gig in (*create expanded go id to gene dictionaries *)
  let extended_go_to_gene_symbol = reverse_dict extended_gsg in
  { fa with
    gene_symbol_to_go = extended_gsg;
    gene_id_to_go = extended_gig;
    go_to_gene_symbol = extended_go_to_gene_symbol;
    go_to_gene_id = extended_go_to_gene_id}


let write_annotations fa path =
  let s_to_go = fa.gene_symbol_to_go in
  let output = Out_channel.create ~append:false path in
  let write_entry k l =
    List.iter l ~f:(fun pkey ->
        let term = Ontology.PKey.get_term fa.ontology pkey in
        Printf.fprintf output "%s\t%s\n" k term.Ontology.Term.id
      )
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
  | Some l -> term_names_of_pkeys fa l
  | None -> []

let create_term_table fa v =
  Ontology.PKey.Table.create fa.ontology v
