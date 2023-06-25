open Core

type t

type domain = Biological_process | Molecular_function | Cellular_component

module Term : sig
  type t = {
    id : string ;
    name : string ;
    namespace : string ;
    is_a : string list ;
  }
  [@@deriving show]

  val of_obo_term :  Obo.term -> t

  val get_parents : t -> string list

  val compare : t -> t -> int
end

val find_term : t -> string -> Term.t option

val expand_term_list : t -> Term.t list -> Term.t list

val expand_id_list : t -> string list -> string list

val filter_terms : t -> string list -> string list

val of_obo : Obo.t -> domain -> (t, string) result

val define_domain : string -> (domain, string) result

val term_names : t -> string String.Map.t
