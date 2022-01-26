type t 

val show : t -> [ `Gene_id_to_go | `Gene_symbol_to_go | `Go_to_gene_symbol | `Go_to_gene_id ] -> string 

val extract_terms : t -> [ `Id of string | `Symbol of string ] -> string list option

val extract_genes : t -> go_id:string -> [ `Id | `Symbol ] -> string list option

val from_gaf_and_ontology : Gaf.t -> Ontology.t -> t   

val propagate_annotations : t -> Ontology.t -> t

val write_annotations : t -> string -> unit


