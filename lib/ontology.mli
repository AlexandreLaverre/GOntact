open Core

type t
type ontology = t

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

(** primary key for terms *)
module PKey : sig
  type t
  type pkey = t
  val compare : t -> t -> int
  val get_term : ontology -> t -> Term.t
  val find_term_by_id : ontology -> string -> t option
  val is_a_transitive_closure : ontology -> t list -> t list

  module Map : Map_intf.S with type Key.t = t

  module Table : sig
    type 'a t
    val create : ontology -> 'a -> 'a t
    val get : 'a t -> pkey -> 'a
    val set : 'a t -> pkey -> 'a -> unit
    val fold : 'a t -> init:'b -> f:(term:Term.t -> value:'a -> 'b -> 'b) -> 'b
  end
end
