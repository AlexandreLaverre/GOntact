open Core
    
type t

val from_ensembl_biomart_file : string -> (t, string) result

val filter_transcript_biotypes : t -> string -> t

val filter_gene_symbols : t -> String.Set.t -> t 

val filter_gene_biotypes : t -> string -> t
  
(*
val show_genes : t -> string

val show_transcripts : t -> string

val filter_by_gene_type : t -> string -> t

val filter_by_transcript_type : t -> string -> t 
*)
