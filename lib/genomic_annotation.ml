open Core

(**************************************************************)

type gene_annot = {
    gene_symbol : string ;
    gene_type : string ;
    chr : string ;
    gene_start : int ;
    gene_end : int ;
    strand : Genomic_interval.strand ;
  }

type exon_annot = {
    transcript_id : string ;
    length : int ;
  }

type transcript_annot = {
    gene_id : string ;
    transcript_type : string ;
    tss_pos : int ;
  }

type t = {
    genes : gene_annot String.Map.t ;
    transcripts : transcript_annot String.Map.t ;
    isoforms : string list String.Map.t ;
    transcript_length : int String.Map.t ;
  }

(**************************************************************)

(* get information from 9th column of GTF *)
(* example: *)
(* 	gene_id "ENSG00000142611"; gene_version "17"; transcript_id "ENST00000511072"; transcript_version "5"; exon_number "1"; gene_name "PRDM16"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; transcript_name "PRDM16-206"; transcript_source "havana"; transcript_biotype "protein_coding"; tag "gencode_basic"; tag "gencode_primary"; transcript_support_level "5"; *)

let info_of_tuple key (k, v) =
  if String.equal k key then Some v else None 

let tuple_of_info_string str =
  let tr = String.strip str in (* remove trailing spaces *)
  let il = String.split tr ~on:'"' in
  match il with
  | [ info_key ; value ; _ ] -> (String.strip info_key, value)
  | _ -> invalid_arg "info should be separated by double quotes (\")"

let tuple_list_of_info_string str =
  let stripped = String.rstrip ~drop:(Char.equal ';') str in
  let info_list = String.split stripped ~on:';' in
  List.filter_map info_list ~f:(function
      | "" -> None
      | s -> Some (tuple_of_info_string s)
    )

(**************************************************************)

(* get transcript annotation from GTF line *)
let transcript_annot_of_gtf_split_line sl =
  match sl with
  | [ _ ; _ ; _ ; transcript_start ; transcript_end ; _ ; transcript_strand ; _ ; info_col ] -> (
    let tss = (
        match transcript_strand with
        | "+" | "1" -> int_of_string transcript_start
        | "-" | "-1" -> int_of_string transcript_end
        | _ -> invalid_arg "wrong strand for gene"
      ) in
    let tuple_list = tuple_list_of_info_string info_col in
    let transcript_id_list = List.filter_map tuple_list ~f:(info_of_tuple "transcript_id") in
    let gene_id_list = List.filter_map tuple_list ~f:(info_of_tuple "gene_id") in
    let transcript_biotype_list = List.filter_map tuple_list ~f:(info_of_tuple "transcript_biotype") in
    match transcript_id_list, gene_id_list, transcript_biotype_list with
    | [ transcript_id ], [ gene_id ], [ transcript_biotype ] ->  (transcript_id,
                                                                  {gene_id = gene_id ; 
                                                                   transcript_type = transcript_biotype ;
                                                                   tss_pos = tss ;
                                                                 })
    | _ -> invalid_arg "couldn't find required info (transcript_id, gene_id, transcript_biotype) in the GTF line"
  )
  | _ -> invalid_arg "GTF line should have 9 tab-separated columns"


(* get exon annotation *)
let exon_annot_of_gtf_split_line sl = 
  match sl with
  | [ _ ; _ ; _ ; exon_start ; exon_end ; _ ; _ ; _ ; info_col ] -> (
    let length = (int_of_string exon_end) - (int_of_string exon_start) +1 in  
    let tuple_list =  tuple_list_of_info_string info_col in
    let transcript_id_list = List.filter_map tuple_list ~f:(info_of_tuple "transcript_id") in
    match transcript_id_list with
    | [ transcript_id ]  -> 
       {transcript_id = transcript_id ; 
        length = length ; }
    | _ -> invalid_arg "couldn't find required info (transcript_id) in the GTF line"
  )
  | _ -> invalid_arg "GTF line should have 9 tab-separated columns"

let gene_annot_of_gtf_split_line sl = 
  match sl with
  | [ chr ; _ ; _ ; gene_start ; gene_end ; _ ; s ; _ ; info_col ] -> (
    let strand = (
        match s with
        | "+" | "1" -> Genomic_interval.Forward
        | "-" | "-1" -> Genomic_interval.Reverse
        | _ -> invalid_arg "wrong strand for gene"
      ) in
    let tuple_list =  tuple_list_of_info_string info_col in
    let gene_id_list = List.filter_map tuple_list ~f:(info_of_tuple "gene_id") in
    let gene_symbol_list = List.filter_map tuple_list ~f:(info_of_tuple "gene_name") in
    let gene_type_list = List.filter_map tuple_list ~f:(info_of_tuple "gene_biotype") in
    match gene_id_list, gene_symbol_list, gene_type_list with
    | [ gene_id ], [ gene_symbol ], [ gene_type ]  -> 
       (gene_id, { gene_symbol = gene_symbol ;
                   gene_type = gene_type ;
                   chr = chr ;
                   gene_start = int_of_string gene_start ;
                   gene_end = int_of_string gene_end ;
                   strand = strand ;})
    | [ gene_id ], [ ], [ gene_type ]  -> 
       (gene_id, { gene_symbol = "NA" ;
                   gene_type = gene_type ;
                   chr = chr ;
                   gene_start = int_of_string gene_start ;
                   gene_end = int_of_string gene_end ;
                   strand = strand ;})
    | _ ->  invalid_argf "couldn't find required info (gene_id, gene_symbol, gene_biotype) in the GTF line %S !" info_col ()
  )
  | _ -> invalid_arg "GTF line should have 9 tab-separated columns"

(**************************************************************)
               
let get_exon_line line =
  let sl = String.split line ~on:'\t' in
  match sl with
  | [ _ ; _ ; "exon" ; _ ; _ ; _ ; _  ; _ ; _ ] -> Some sl
  | _ -> None


 let get_transcript_line line =
  let sl = String.split line ~on:'\t' in
  match sl with
  | [ _ ; _ ; "transcript" ; _ ; _ ; _ ; _  ; _ ; _ ] -> Some sl
  |  _ -> None


let get_gene_line line =
  let sl = String.split line ~on:'\t' in
  match sl with
  | [ _ ; _ ; "gene" ; _ ; _ ; _ ; _  ; _ ; _ ] -> Some sl
  | _ ->  None



let isoforms_of_transcripts tx =
  Map.map tx ~f:(fun t -> t.gene_id) (*simplify dictionary - tx id - gene id *)
  |> Map.to_alist (*transform to list of tuples*)
  |> List.map ~f:(fun (x, y) -> (y, x)) (*reverse tuple ordre to create dictionary *)
  |> String.Map.of_alist_multi  (*dictionary has gene ids as keys, values are list of transcript ids *)

let of_gtf_file path =
  let open Let_syntax.Result in
  let* lines = Utils.read_lines path in
  print_endline "piu0" ;
  match lines with
  | _ :: _ -> (
    let exon_lines = List.filter_map lines ~f:get_exon_line in
    print_endline "piu1" ;
    let exon_list = List.map exon_lines ~f:exon_annot_of_gtf_split_line in
    print_endline "piu2" ;
    let exon_map =  List.fold exon_list ~init:String.Map.empty ~f:(fun acc r ->
                        Map.add_multi acc ~key:r.transcript_id ~data:r
                      ) in
    print_endline "piu3" ;
    let transcript_length =  Map.map exon_map ~f:(fun exonlist -> List.fold exonlist ~init:0 ~f:(fun x y -> x + y.length)) in 
    print_endline "piu4" ;
    let transcript_lines = List.filter_map lines ~f:get_transcript_line in
    print_endline "piu5" ;
    let transcript_list = List.map transcript_lines ~f:transcript_annot_of_gtf_split_line in
    print_endline "piu6" ;
    let transcripts = String.Map.of_alist_exn transcript_list in
    print_endline "piu7" ;
    let gene_lines = List.filter_map lines ~f:get_gene_line in
    print_endline "piu8" ;
    let gene_list = List.map gene_lines ~f:gene_annot_of_gtf_split_line in
    print_endline "piu9" ;
    let genes = String.Map.of_alist_exn gene_list in
    print_endline "piu10" ;
    let isoforms = isoforms_of_transcripts transcripts in
    print_endline "piu11" ;
    Ok {genes ; transcripts; isoforms ; transcript_length}
  )
  | [] -> Error "File is empty"

(**************************************************************)

let gene_symbol ga id =
  match Map.find ga.genes id with
  | Some g -> Some g.gene_symbol
  | None -> None

let gene_symbol_exn ga id =
  let g = Option.value_exn (Map.find ga.genes id) in
  g.gene_symbol

let filter_transcript_biotypes ga biotype =
  let genes = ga.genes in
  let transcripts = ga.transcripts in
  let transcript_length = ga.transcript_length in
  let filtered_transcripts = Map.filter transcripts ~f:(fun x -> String.equal x.transcript_type biotype) in (*select transcripts with the good biotype*)
  let filtered_gene_set = Map.fold filtered_transcripts ~init:String.Set.empty ~f:(fun ~key:_ ~data:tx acc -> Set.add acc tx.gene_id) in
  let filtered_genes = Map.filter_keys genes ~f:(Set.mem filtered_gene_set) in
  let filtered_isoforms = isoforms_of_transcripts filtered_transcripts in
  {genes = filtered_genes ; transcripts = filtered_transcripts ; isoforms = filtered_isoforms ; transcript_length}

let filter_gene_biotypes ga biotype =
  let genes = ga.genes in
  let transcripts = ga.transcripts in
  let filtered_genes = Map.filter genes ~f:(fun x -> String.equal x.gene_type biotype) in (*select genes with the good biotype*)
  let filtered_gene_set = String.Set.of_list (Map.keys filtered_genes) in
  let filtered_transcripts = Map.filter transcripts ~f:(fun x -> Set.mem filtered_gene_set x.gene_id) in
  let filtered_isoforms = isoforms_of_transcripts filtered_transcripts in
  {genes = filtered_genes ; transcripts = filtered_transcripts ; isoforms = filtered_isoforms ; transcript_length = ga.transcript_length}

let filter_gene_symbols ga symbol_set =
  (*symbol set is a set of strings = gene symbols*)
  let filtered_genes = Map.filter ga.genes ~f:(fun g -> Set.mem symbol_set g.gene_symbol) in (*we keep genes that have these symbols*)
  let gene_set = String.Set.of_list (Map.keys filtered_genes) in
  let filtered_transcripts = Map.filter ga.transcripts ~f:(fun tx -> Set.mem gene_set tx.gene_id) in
  let filtered_isoforms = Map.filter_keys ga.isoforms ~f:(fun gid -> Set.mem gene_set gid) in
  {genes = filtered_genes ; transcripts = filtered_transcripts ; isoforms = filtered_isoforms ; transcript_length = ga.transcript_length}

let filter_chromosomes ga chr_set =
  let filtered_genes = Map.filter ga.genes ~f:(fun g -> Set.mem chr_set g.chr) in (*we keep genes that are on these chromosomes*)
  let gene_set = String.Set.of_list (Map.keys filtered_genes) in
  let filtered_transcripts = Map.filter ga.transcripts ~f:(fun tx -> Set.mem gene_set tx.gene_id) in
  let filtered_isoforms = Map.filter_keys ga.isoforms ~f:(fun gid -> Set.mem gene_set gid) in
  {genes = filtered_genes ; transcripts = filtered_transcripts ; isoforms = filtered_isoforms  ; transcript_length = ga.transcript_length}

let remove_duplicated_gene_symbols ga =
  let tuples = Map.fold ga.genes ~init:[]  ~f:(fun ~key:gene_id ~data:gene acc -> ((gene.gene_symbol, gene_id) :: acc)) in
  let symbol_id = String.Map.of_alist_multi tuples in
  let unique_genes = Map.fold symbol_id ~init:[] ~f:(fun ~key:_ ~data acc -> (
        match data with
        | [ gid ] -> gid :: acc  (*keep only those genes that are unique*)
        | _ -> acc
      )
    ) in
  let filtered_gene_tuples = List.map unique_genes ~f:(fun g -> (g, (Map.find_exn ga.genes g))) in
  let filtered_genes = String.Map.of_alist_exn filtered_gene_tuples in
  let gene_set = String.Set.of_list (Map.keys filtered_genes) in
  let filtered_transcripts = Map.filter ga.transcripts ~f:(fun tx -> Set.mem gene_set tx.gene_id) in
  let filtered_isoforms = Map.filter_keys ga.isoforms ~f:(fun gid -> Set.mem gene_set gid) in
  {genes = filtered_genes ; transcripts = filtered_transcripts ; isoforms = filtered_isoforms ; transcript_length = ga.transcript_length}


let compare_isoforms txlength t1 t2 =
  (*two isoforms from the same gene *)
  let len1 = Map.find_exn txlength t1 in (*t1 and t2 are necessarily in the transcript dictionary*)
  let len2 = Map.find_exn txlength t2 in
  compare len1 len2

let identify_major_isoforms ga =
  let gene_list = Map.keys ga.genes in
  let isoforms = ga.isoforms in
  let transcript_length = ga.transcript_length in
  let find_major_isoform isolist =
    Option.value_exn (List.max_elt isolist ~compare:(compare_isoforms transcript_length))
  in
  let major_list = List.map gene_list ~f:(fun g -> (g, find_major_isoform (Map.find_exn isoforms g))) in
  String.Map.of_alist_exn major_list

let identify_major_isoforms_symbols ga =
  let gene_list = Map.keys ga.genes in
  let isoforms = ga.isoforms in
  let transcript_length = ga.transcript_length in
  let find_major_isoform isolist =
    Option.value_exn (List.max_elt isolist ~compare:(compare_isoforms transcript_length))
  in
  let major_list = List.map gene_list ~f:(fun g -> ((gene_symbol_exn ga g), find_major_isoform (Map.find_exn isoforms g))) in
  String.Map.of_alist_exn major_list

let write_major_isoforms mi path =
  Out_channel.with_file path ~append:false ~f:(fun output ->
      Printf.fprintf output "GeneID\tTranscript\n" ;
      Map.iteri mi  ~f:(fun ~key ~data -> Printf.fprintf output "%s\t%s\n" key data))

let major_isoform_tss ga ~major_isoforms:iso =
  let ginfo = ga.genes in
  let txinfo = ga.transcripts in
  let tss_pos = Map.map iso ~f:(fun tx -> (let ti = Map.find_exn txinfo tx in ti.tss_pos)) in
  let interval_of_gene g =
    let gi = Map.find_exn ginfo g in
    let tss = Map.find_exn tss_pos g in
    Genomic_interval.make ~id:g gi.chr tss tss gi.strand
  in
  let int_list = List.map (Map.keys iso) ~f:interval_of_gene in
  Genomic_interval_collection.of_interval_list int_list

let all_tss_intervals ga extend =
  let ginfo = ga.genes in
  let txinfo = ga.transcripts in
  let tss_pos = Map.map txinfo ~f:(fun tx -> tx.tss_pos) in
  let interval_of_transcript tx  =
    let ti = Map.find_exn txinfo tx in
    let g = ti.gene_id in
    let gi = Map.find_exn ginfo g in
    let gs = gi.gene_symbol in
    let new_id = Printf.sprintf "%s:%s" g gs in
    let tss = Map.find_exn tss_pos tx in
    let max_extend = max extend 1 in
    Genomic_interval.make ~id:new_id gi.chr (tss - max_extend) (tss + max_extend) gi.strand
  in
  let int_list = List.map (Map.keys txinfo) ~f:interval_of_transcript in
  Genomic_interval_collection.of_interval_list int_list

let compute_cis_distances element_genes ~element_map:em ~gene_annotation:ga ~major_isoforms:mi =
  (* we assume that pairs of elements-genes are always in cis *)
  let txinfo = ga.transcripts in
  let tss_pos = Map.map txinfo ~f:(fun tx -> tx.tss_pos) in
  let compute_one_distance element_id gene_symbol =
    let txid = Map.find_exn mi gene_symbol in (*identifier of the major isoform*)
    let this_tss_pos = float_of_int (Map.find_exn tss_pos txid) in
    let gint = Map.find_exn em element_id in
    let start_pos = float_of_int (Genomic_interval.start_pos gint) in
    let end_pos = float_of_int (Genomic_interval.end_pos gint) in
    let mid_pos = (start_pos +. end_pos) /. 2. in
    let dist = Float.abs (mid_pos -. this_tss_pos) in
    (gene_symbol, dist)
  in
  Map.mapi element_genes ~f:(fun ~key ~data -> (List.map data ~f:(fun g -> compute_one_distance key g)))

let write_distance_elements ~dist_elements path =
  Out_channel.with_file path ~append:false ~f:(fun output ->
      Printf.fprintf output "ElementID\tGeneSymbol\tDistance\n" ;
      Map.iteri dist_elements  ~f:(fun ~key ~data -> (List.iter data ~f:(fun (gs, dist) -> Printf.fprintf output "%s\t%s\t%f\n" key gs dist))))

let number_of_genes ga =
  Map.length ga.genes

(**************************************************************)

(* get chromosome size from gene coordinate list *)
(*let chromosome_sizes_of_annot ga extend =
  let gene_map = ga.genes in
  let gene_ends_map =  Map.fold gene_map ~init:String.Map.empty ~f:(fun acc r ->
                        Map.add_multi acc ~key:r.chr ~data:r.gene_end
*)
