type t
val of_sorted_lists_unsafe : Ontology.PKey.t list list -> t
val to_sorted_list : t -> Ontology.PKey.t list
val iter : t -> f:(Ontology.PKey.t -> unit) -> unit
val union : t -> t -> t
