type t
[@@deriving show]

type interval_comparison = Smaller_chr | Larger_chr | Smaller_no_overlap | Smaller_overlap | Equal | Larger_overlap | Larger_no_overlap

type strand = Forward | Reverse | Unstranded
[@@deriving show]

val string_of_strand : strand -> string
  
val chr : t -> string

val id : t -> string 

val start_pos : t -> int

val end_pos : t -> int 

val strand : t -> strand
  
val make : ?id:string -> string -> int -> int -> strand ->  t

val check_overlap : t -> t -> interval_comparison

val compare_intervals : t -> t -> int 

val merge : t -> t -> t option
    
val intersect : t -> t -> [ `Stranded | `Unstranded ]-> bool
