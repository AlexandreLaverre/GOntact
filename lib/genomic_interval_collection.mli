type t
[@@deriving show]

val from_bed_file : string -> t

val merge_coordinates : t -> t
(*val sort_by_coordinate : t -> t

  val extract_intersection : t -> t -> t 
*)
