type t
[@@deriving show]

val from_ensembl_biomart_file : string -> (t, string) result

val show_genes : t -> string

val show_transcripts : t -> string
(*
val filter_by_gene_type : t -> string -> t

val filter_by_transcript_type : t -> string -> t 
*)
