type t
[@@deriving show]

val chr : t -> string

val id : t -> string 

val start_pos : t -> int

val end_pos : t -> int 

val make : ?id:string -> string -> int -> int -> t

val compare : t -> t -> int

val merge : t -> t -> t option
    
val intersect : t -> t -> bool
