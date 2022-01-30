(** GAF format specification: http://geneontology.org/docs/go-annotation-file-gaf-format-2.1/ *)

type gene_annotation = {
  gene_id : string ;
  gene_symbol : string ;
  qualifier : string ;
  go_id : string ;
}
[@@deriving show]

type t = gene_annotation list
[@@deriving show]

val of_gaf_file : string -> (t, string) result


