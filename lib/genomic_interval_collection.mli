open Core

type t
[@@deriving show]

type format = Base1 | Base0 

val of_bed_file : string -> strip_chr:bool -> format:format -> t

val of_chr_size_file : string -> strip_chr:bool -> t

val of_interval_list : Genomic_interval.t list -> t

val interval_list : t -> Genomic_interval.t list 

val chr_set :  t -> String.Set.t 

val merge_coordinates : t -> t

val intersect : t -> t -> string list String.Map.t 
 
val map : t -> f:(Genomic_interval.t -> Genomic_interval.t) -> t

val iter : t -> f:(Genomic_interval.t -> unit) -> unit

val reverse_sort_by_coordinate : t -> t

val write_output : t -> string -> append:bool  -> unit 

