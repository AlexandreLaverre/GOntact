type term = {
  id : string ;
  name : string ;
  namespace : string ;
  is_a : string list ;
}

type t = term list
    
(*val from_file : string -> t*)
  
