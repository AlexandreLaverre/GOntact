open Core
    
type t

val of_ensembl_biomart_file : string -> (t, string) result

val filter_transcript_biotypes : t -> string -> t

val filter_gene_symbols : t -> String.Set.t -> t 

val filter_gene_biotypes : t -> string -> t

val filter_chromosomes : t -> String.Set.t -> t

val identify_major_isoforms : t -> string String.Map.t

val major_isoform_tss : t -> major_isoforms: string String.Map.t -> Genomic_interval_collection.t 

val all_tss_intervals : t -> int -> Genomic_interval_collection.t

val gene_symbol : t -> string -> string option

val gene_symbol_exn : t -> string -> string 

(*
val show_genes : t -> string

val show_transcripts : t -> string

val filter_by_gene_type : t -> string -> t

val filter_by_transcript_type : t -> string -> t 
*)
