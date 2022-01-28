open Core

type gene_annot = {
  gene_symbol : string ;
  gene_type : string ;
  chr : string ;
  strand : string ; 
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

let gene_annot_from_line line header =
  let sl = String.split line ~on:'\t' in
  let al = Array.of_list sl in
  let gene_id = al.(header.gene_id_index) in
  (gene_id , 
  {gene_symbol = al.(header.gene_symbol_index) ;
   gene_type = al.(header.gene_type_index) ;
   chr = al.(header.chr_index) ;
   strand = al.(header.strand_index)})
  
let transcript_annot_from_line line header =
  let sl = String.split line ~on:'\t' in
  let al = Array.of_list sl in
  let transcript_id = al.(header.transcript_id_index) in
  (transcript_id, 
   {gene_id = al.(header.gene_id_index) ;
    transcript_type = al.(header.transcript_type_index) ; 
    appris_class = al.(header.appris_index) ;
    tss_pos = int_of_string (al.(header.tss_pos_index)) ;
    length = int_of_string (al.(header.transcript_length_index))})

let isoforms_from_transcripts tx =
  String.Map.map tx ~f:(fun t -> t.gene_id) (*simplify dictionary - tx id - gene id *)  
  |> String.Map.to_alist (*transform to list of tuples*)
  |> List.map ~f:(fun (x, y) -> (y, x)) (*reverse tuple ordre to create dictionary *)
  |> String.Map.of_alist_multi  (*dictionary has gene ids as keys, values are list of transcript ids *)
                                              
let compare_gene_tuples g1 g2 =
  let (g1id, _) = g1 in
  let (g2id, _) = g2 in
  String.compare g1id g2id
    
let from_ensembl_biomart_file path =
  let lines = In_channel.read_lines path in
  match lines with
  | h :: t -> (
      let open Let_syntax.Result in
      let+ header = extract_ensembl_biomart_header h in
      let gene_list = List.map t ~f:(fun line -> gene_annot_from_line line header) in
      let dedup_gene_list = List.dedup_and_sort ~compare:compare_gene_tuples gene_list in  
      let transcript_list = List.map t ~f:(fun line -> transcript_annot_from_line line header) in
      let genes = String.Map.of_alist_exn dedup_gene_list in
      let transcripts = String.Map.of_alist_exn transcript_list in
      let isoforms = isoforms_from_transcripts transcripts in
      {genes ; transcripts; isoforms}
    )
  | [] -> Error "File is empty."


let filter_transcript_biotypes ga biotype =
  let genes = ga.genes in 
  let transcripts = ga.transcripts in
  let filtered_transcripts = String.Map.filter transcripts ~f:(fun x -> String.equal x.transcript_type biotype) in (*select transcripts with the good biotype*)
  let filtered_gene_set = String.Map.fold filtered_transcripts ~init:String.Set.empty ~f:(fun ~key:_ ~data:tx acc -> String.Set.add acc tx.gene_id) in  
  let filtered_genes = String.Map.filter_keys genes ~f:(String.Set.mem filtered_gene_set) in
  let filtered_isoforms = isoforms_from_transcripts filtered_transcripts in
  {genes = filtered_genes ; transcripts = filtered_transcripts ; isoforms = filtered_isoforms}

let filter_gene_biotypes ga biotype =
  let genes = ga.genes in 
  let transcripts = ga.transcripts in
  let filtered_genes = String.Map.filter genes ~f:(fun x -> String.equal x.gene_type biotype) in (*select genes with the good biotype*)
  let filtered_gene_set = String.Set.of_list (String.Map.keys filtered_genes) in
  let filtered_transcripts = String.Map.filter transcripts ~f:(fun x -> String.Set.mem filtered_gene_set x.gene_id) in
  let filtered_isoforms = isoforms_from_transcripts filtered_transcripts in
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

(*
let show_genes annot =
  [%show: gene_annot list] annot.genes

let show_transcripts annot =
  [%show: transcript_annot list] annot.transcripts
  *)
