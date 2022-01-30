type t
[@@deriving show]

val chr : t -> string

val id : t -> string 

val start_pos : t -> int

val end_pos : t -> int 

val strand : t -> string
  
val make : ?id:string -> string -> int -> int -> string ->  t

val compare : t -> t -> int

val merge : t -> t -> t option
    
val intersect : t -> t -> [ `Stranded | `Unstranded ]-> bool
