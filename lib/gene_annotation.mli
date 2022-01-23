type t 

val extract_terms : t -> [ `Id of string | `Symbol of string ] -> string list option

val extract_genes : t -> go_id:string -> [ `Id | `Symbol ] -> string list option

val from_gaf : Gaf.t -> t   
