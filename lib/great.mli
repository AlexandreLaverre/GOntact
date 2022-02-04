type regulatory_domain 

val genomic_interval_collection : regulatory_domain list -> Genomic_interval_collection.t 

val basal_plus_extension_domains_one_chr : chr:string -> chromosome_size:int -> genomic_annotation:Genomic_annotation.t ->  upstream:int -> downstream:int -> extend:int -> regulatory_domain list

val basal_plus_extension_domains : chromosome_sizes:Genomic_interval_collection.t -> genomic_annotation:Genomic_annotation.t ->  upstream:int -> downstream:int -> extend:int -> regulatory_domain list

(*
val go_term_frequencies : t ->
*)
