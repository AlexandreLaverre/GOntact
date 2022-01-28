open Core

type t
[@@deriving show]

val from_bed_file : string -> t

val from_chr_size_file : string -> t

val chr_set :  t -> String.Set.t 

val merge_coordinates : t -> t

(*
  val extract_intersection : t -> t -> t 
*)
