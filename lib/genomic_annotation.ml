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
      {genes ; transcripts}
    )
  | [] -> Error "File is empty."


let filter_transcript_biotypes ga biotype =
  let genes = ga.genes in 
  let transcripts = ga.transcripts in
  let filtered_transcripts = String.Map.filter transcripts ~f:(fun x -> String.equal x.transcript_type biotype) in
  let tx_info_list = String.Map.data filtered_transcripts in
  let gene_id_list = List.map tx_info_list ~f:(fun x -> x.gene_id) in 
  let gene_tuple_list = List.map gene_id_list ~f:(fun g -> (g, String.Map.find_exn genes g)) in
  let filtered_genes = String.Map.of_alist_exn gene_tuple_list in
  {genes = filtered_genes ; transcripts = filtered_transcripts}

(*
let show_genes annot =
  [%show: gene_annot list] annot.genes

let show_transcripts annot =
  [%show: transcript_annot list] annot.transcripts
  *)
