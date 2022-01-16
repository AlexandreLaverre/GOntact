type term = {
  id : string ;
  name : string ;
  namespace : string ;
  is_a : string list ;
}
[@@deriving show]

type t = term list
[@@deriving show]

val from_file : string -> (t, string) result
