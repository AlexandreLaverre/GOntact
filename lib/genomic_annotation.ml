open Core

type gene_annot = {
  gene_symbol : string ;
  gene_type : string ;
  chr : string ;
  strand : Genomic_interval.strand ; 
}

type transcript_annot = {
  gene_id : string ;
  transcript_type : string ; 
  appris_class : string ;
  tss_pos : int ;
  length : int ;
}

type t = {
  genes : gene_annot String.Map.t ;
  transcripts : transcript_annot String.Map.t ;
  isoforms : string list String.Map.t ;
}

type biomart_header = {
  gene_id_index : int ;
  gene_type_index : int ;
  gene_symbol_index : int ;
  transcript_id_index : int ;
  transcript_type_index : int ;
  appris_index : int ;
  chr_index : int ;
  strand_index : int ;
  tss_pos_index : int ;
  transcript_length_index : int ;
}

let gene_symbol ga id =
  match String.Map.find ga.genes id with
  | Some g -> Some g.gene_symbol
  | None -> None

let gene_symbol_exn ga id =
  let g = Option.value_exn (String.Map.find ga.genes id) in
  g.gene_symbol

let extract_ensembl_biomart_header h =
  let sh = String.split h ~on:'\t' in
  let index str =
    match List.findi sh ~f:(fun _ elt -> String.equal elt str) with
    | Some (i, _) -> Ok i
    | None -> Error (Printf.sprintf "Cannot find column %s." str)
  in
  let header_list = ["Gene stable ID" ; "Gene type" ; "Gene name" ; "Transcript stable ID" ; "Transcript type" ; "APPRIS annotation" ; "Chromosome/scaffold name" ; "Strand" ; "Transcription start site (TSS)" ; "Transcript length (including UTRs and CDS)" ] in
  let index_res = List.map header_list ~f:(fun x -> index x) in
  let indexes = Result.all index_res in
  match indexes with
  | Ok [ gene_id_index ; gene_type_index ; gene_symbol_index ; transcript_id_index ; transcript_type_index  ; appris_index ; chr_index ; strand_index ;  tss_pos_index ; transcript_length_index ] ->  Ok { gene_id_index ; gene_type_index ; gene_symbol_index ; transcript_id_index ; transcript_type_index  ; appris_index ; chr_index ; strand_index ;  tss_pos_index ; transcript_length_index }  
  | _ -> Error "Ensembl BioMart file does not have all necessary columns."

let gene_annot_of_line line header =
  let sl = String.split line ~on:'\t' in
  let al = Array.of_list sl in
  let s = al.(header.strand_index) in
  let strand = (
    match s with
    | "+" | "1" -> Genomic_interval.Forward
    | "-" | "-1" -> Genomic_interval.Reverse
    | _ -> invalid_arg "wrong strand for gene"
  ) in  
  let gene_id = al.(header.gene_id_index) in
  (gene_id, 
  {gene_symbol = al.(header.gene_symbol_index) ;
   gene_type = al.(header.gene_type_index) ;
   chr = al.(header.chr_index) ;
   strand})
  
let transcript_annot_of_line line header =
  let sl = String.split line ~on:'\t' in
  let al = Array.of_list sl in
  let transcript_id = al.(header.transcript_id_index) in
  (transcript_id, 
   {gene_id = al.(header.gene_id_index) ;
    transcript_type = al.(header.transcript_type_index) ; 
    appris_class = al.(header.appris_index) ;
    tss_pos = int_of_string (al.(header.tss_pos_index)) ;
    length = int_of_string (al.(header.transcript_length_index))})

let isoforms_of_transcripts tx =
  String.Map.map tx ~f:(fun t -> t.gene_id) (*simplify dictionary - tx id - gene id *)  
  |> String.Map.to_alist (*transform to list of tuples*)
  |> List.map ~f:(fun (x, y) -> (y, x)) (*reverse tuple ordre to create dictionary *)
  |> String.Map.of_alist_multi  (*dictionary has gene ids as keys, values are list of transcript ids *)
                                              
let compare_gene_tuples g1 g2 =
  let (g1id, _) = g1 in
  let (g2id, _) = g2 in
  String.compare g1id g2id
    
let of_ensembl_biomart_file path =
  let lines = In_channel.read_lines path in
  match lines with
  | h :: t -> (
      let open Let_syntax.Result in
      let+ header = extract_ensembl_biomart_header h in
      let gene_list = List.map t ~f:(fun line -> gene_annot_of_line line header) in
      let dedup_gene_list = List.dedup_and_sort ~compare:compare_gene_tuples gene_list in  
      let transcript_list = List.map t ~f:(fun line -> transcript_annot_of_line line header) in
      let genes = String.Map.of_alist_exn dedup_gene_list in
      let transcripts = String.Map.of_alist_exn transcript_list in
      let isoforms = isoforms_of_transcripts transcripts in
      {genes ; transcripts; isoforms}
    )
  | [] -> Error "File is empty."


let filter_transcript_biotypes ga biotype =
  let genes = ga.genes in 
  let transcripts = ga.transcripts in
  let filtered_transcripts = String.Map.filter transcripts ~f:(fun x -> String.equal x.transcript_type biotype) in (*select transcripts with the good biotype*)
  let filtered_gene_set = String.Map.fold filtered_transcripts ~init:String.Set.empty ~f:(fun ~key:_ ~data:tx acc -> String.Set.add acc tx.gene_id) in  
  let filtered_genes = String.Map.filter_keys genes ~f:(String.Set.mem filtered_gene_set) in
  let filtered_isoforms = isoforms_of_transcripts filtered_transcripts in
  {genes = filtered_genes ; transcripts = filtered_transcripts ; isoforms = filtered_isoforms}

let filter_gene_biotypes ga biotype =
  let genes = ga.genes in 
  let transcripts = ga.transcripts in
  let filtered_genes = String.Map.filter genes ~f:(fun x -> String.equal x.gene_type biotype) in (*select genes with the good biotype*)
  let filtered_gene_set = String.Set.of_list (String.Map.keys filtered_genes) in
  let filtered_transcripts = String.Map.filter transcripts ~f:(fun x -> String.Set.mem filtered_gene_set x.gene_id) in
  let filtered_isoforms = isoforms_of_transcripts filtered_transcripts in
  {genes = filtered_genes ; transcripts = filtered_transcripts ; isoforms = filtered_isoforms}
  
let filter_gene_symbols ga symbol_set =
  (*symbol set is a set of strings = gene symbols*)
  let filtered_genes = String.Map.filter ga.genes ~f:(fun g -> String.Set.mem symbol_set g.gene_symbol) in (*we keep genes that have these symbols*)
  let gene_set = String.Set.of_list (String.Map.keys filtered_genes) in
  let filtered_transcripts = String.Map.filter ga.transcripts ~f:(fun tx -> String.Set.mem gene_set tx.gene_id) in
  let filtered_isoforms = String.Map.filter_keys ga.isoforms ~f:(fun gid -> String.Set.mem gene_set gid) in
  {genes = filtered_genes ; transcripts = filtered_transcripts ; isoforms = filtered_isoforms}

let filter_chromosomes ga chr_set =
  let filtered_genes = String.Map.filter ga.genes ~f:(fun g -> String.Set.mem chr_set g.chr) in (*we keep genes that are on these chromosomes*)
  let gene_set = String.Set.of_list (String.Map.keys filtered_genes) in
  let filtered_transcripts = String.Map.filter ga.transcripts ~f:(fun tx -> String.Set.mem gene_set tx.gene_id) in
  let filtered_isoforms = String.Map.filter_keys ga.isoforms ~f:(fun gid -> String.Set.mem gene_set gid) in
  {genes = filtered_genes ; transcripts = filtered_transcripts ; isoforms = filtered_isoforms}

let remove_duplicated_gene_symbols ga =
  let tuples = String.Map.fold ga.genes ~init:[]  ~f:(fun ~key:gene_id ~data:gene acc -> ((gene.gene_symbol, gene_id) :: acc)) in
  let symbol_id = String.Map.of_alist_multi tuples in
  let unique_genes = String.Map.fold symbol_id ~init:[] ~f:(fun ~key:_ ~data acc -> (
        match data with
        | [ gid ] -> gid :: acc  (*keep only those genes that are unique*)
        | _ -> acc
      ) 
    ) in
  let filtered_gene_tuples = List.map unique_genes ~f:(fun g -> (g, (String.Map.find_exn ga.genes g))) in
  let filtered_genes = String.Map.of_alist_exn filtered_gene_tuples in
  let gene_set = String.Set.of_list (String.Map.keys filtered_genes) in
  let filtered_transcripts = String.Map.filter ga.transcripts ~f:(fun tx -> String.Set.mem gene_set tx.gene_id) in
  let filtered_isoforms = String.Map.filter_keys ga.isoforms ~f:(fun gid -> String.Set.mem gene_set gid) in
  {genes = filtered_genes ; transcripts = filtered_transcripts ; isoforms = filtered_isoforms}
            
      
let compare_isoforms txinfo t1 t2 =
  (*two isoforms from the same gene *)
  let info1 = String.Map.find_exn txinfo t1 in (*t1 and t2 are necessarily in the transcript dictionary*)
  let info2 = String.Map.find_exn txinfo t2 in
  let appris1 = info1.appris_class in
  let appris2 = info2.appris_class in
  if String.equal appris1 appris2 then (
    let (len1, len2) = (info1.length, info2.length) in
    compare len1 len2
  )
  else (
    match (appris1, appris2) with
    | ("principal1", _) -> 1
    | (_, "principal1") -> -1
    | ("principal2", _) -> 1
    | (_, "principal2") -> -1
    | ("principal3", _) -> 1
    | (_, "principal3") -> -1
    | ("principal4", _) -> 1
    | (_, "principal4") -> -1
    | ("principal5", _) -> 1
    | (_, "principal5") -> -1
    | ("alternative1", _) -> 1
    | (_, "alternative1") -> -1
    | ("alternative2", _) -> 1
    | (_, "alternative2") -> -1
    | _ -> 0
  )

let identify_major_isoforms ga =
  let gene_list = String.Map.keys ga.genes in
  let isoforms = ga.isoforms in
  let transcripts = ga.transcripts in
  let find_major_isoform isolist =
    Option.value_exn (List.max_elt isolist ~compare:(compare_isoforms transcripts))
  in
  let major_list = List.map gene_list ~f:(fun g -> (g, find_major_isoform (String.Map.find_exn isoforms g))) in
  String.Map.of_alist_exn major_list

let identify_major_isoforms_symbols ga =
  let gene_list = String.Map.keys ga.genes in
  let isoforms = ga.isoforms in
  let transcripts = ga.transcripts in
  let find_major_isoform isolist =
    Option.value_exn (List.max_elt isolist ~compare:(compare_isoforms transcripts))
  in
  let major_list = List.map gene_list ~f:(fun g -> ((gene_symbol_exn ga g), find_major_isoform (String.Map.find_exn isoforms g))) in
  String.Map.of_alist_exn major_list

let write_major_isoforms mi path =
  Out_channel.with_file path ~append:false ~f:(fun output ->
      Printf.fprintf output "GeneID\tTranscript\n" ;
      String.Map.iteri mi  ~f:(fun ~key ~data -> Printf.fprintf output "%s\t%s\n" key data))

let major_isoform_tss ga ~major_isoforms:iso =
  let ginfo = ga.genes in 
  let txinfo = ga.transcripts in
  let tss_pos = String.Map.map iso ~f:(fun tx -> (let ti = String.Map.find_exn txinfo tx in ti.tss_pos)) in  
  let interval_of_gene g =
    let gi = String.Map.find_exn ginfo g in
    let tss = String.Map.find_exn tss_pos g in
    Genomic_interval.make ~id:g gi.chr tss tss gi.strand
  in
  let int_list = List.map (String.Map.keys iso) ~f:interval_of_gene in
  Genomic_interval_collection.of_interval_list int_list

let all_tss_intervals ga extend = 
  let ginfo = ga.genes in 
  let txinfo = ga.transcripts in
  let tss_pos = String.Map.map txinfo ~f:(fun tx -> tx.tss_pos) in  
  let interval_of_transcript tx  =
    let ti = String.Map.find_exn txinfo tx in
    let g = ti.gene_id in 
    let gi = String.Map.find_exn ginfo g in
    let gs = gi.gene_symbol in
    let new_id = Printf.sprintf "%s:%s" g gs in
    let tss = String.Map.find_exn tss_pos tx in
    let max_extend = max extend 1 in 
    Genomic_interval.make ~id:new_id gi.chr (tss - max_extend) (tss + max_extend) gi.strand
  in
  let int_list = List.map (String.Map.keys txinfo) ~f:interval_of_transcript in
  Genomic_interval_collection.of_interval_list int_list
  
let compute_cis_distances element_genes ~element_map:em ~gene_annotation:ga ~major_isoforms:mi =
  (* we assume that pairs of elements-genes are always in cis *) 
  let txinfo = ga.transcripts in
  let tss_pos = String.Map.map txinfo ~f:(fun tx -> tx.tss_pos) in
  let compute_one_distance element_id gene_symbol =
    let txid = String.Map.find_exn mi gene_symbol in (*identifier of the major isoform*)
    let this_tss_pos = float_of_int (String.Map.find_exn tss_pos txid) in
    let gint = String.Map.find_exn em element_id in 
    let start_pos = float_of_int (Genomic_interval.start_pos gint) in
    let end_pos = float_of_int (Genomic_interval.end_pos gint) in
    let mid_pos = (start_pos +. end_pos) /. 2. in
    let dist = Float.abs (mid_pos -. this_tss_pos) in
    (gene_symbol, dist)
  in
  String.Map.mapi element_genes ~f:(fun ~key ~data -> (List.map data ~f:(fun g -> compute_one_distance key g)))

let write_distance_elements ~dist_elements path =
  Out_channel.with_file path ~append:false ~f:(fun output ->
      Printf.fprintf output "ElementID\tGeneSymbol\tDistance\n" ;
      String.Map.iteri dist_elements  ~f:(fun ~key ~data -> (List.iter data ~f:(fun (gs, dist) -> Printf.fprintf output "%s\t%s\t%f\n" key gs dist))))

let number_of_genes ga =
  String.Map.length ga.genes




