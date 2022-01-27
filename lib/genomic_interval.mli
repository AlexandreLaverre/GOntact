type t
[@@deriving show]
  
val make : ?id:string -> string -> int -> int -> t

val compare : t -> t -> int

val merge : t -> t -> t option
    
val intersect : t -> t -> bool
