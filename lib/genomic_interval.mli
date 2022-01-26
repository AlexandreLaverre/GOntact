type t

val make : ?id:string -> string -> int -> int -> t
  
val intersect : t -> t -> bool
