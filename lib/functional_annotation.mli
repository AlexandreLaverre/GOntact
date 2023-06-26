(** Two-way relation between a set of gene identifiers and GO terms *)

type t

(* val show : t -> [ `Gene_id_to_go | `Gene_symbol_to_go | `Go_to_gene_symbol | `Go_to_gene_id ] -> string *)

val extract_terms : t -> [ `Id of string | `Symbol of string ] -> Ontology.PKey.t list option

val extract_terms_exn : t -> [ `Id of string | `Symbol of string ] -> Ontology.PKey.t list

val term_names_of_pkeys : t -> Ontology.PKey.t list -> string list

val extract_genes : t -> go_id:string -> [ `Id | `Symbol ] -> string list option

val of_gaf_and_ontology : Gaf.t -> Ontology.t -> t

val propagate_annotations : t -> Ontology.t -> t

val write_annotations : t -> string -> unit

val gene_symbols : t -> Core.String.Set.t

val go_list_of_gene_symbol : t -> string -> string list

val create_term_table : t -> 'a -> 'a Ontology.PKey.Table.t
