type term = {
  id : string ;
  name : string ;
  namespace : string ;
  is_a : string list ;
}
[@@deriving show]

type t = term list
[@@deriving show]

val of_obo_file : string -> (t, string) result

val filter_namespace : t-> namespace:string -> t
