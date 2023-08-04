open Core

type t = {
  bait_chr : string ;
  bait_start : int ;
  bait_end : int ;
  bait_name : string ;
  otherEnd_chr : string ;
  otherEnd_start : int ;
  otherEnd_end : int ;
  otherEnd_name : string ;
  n_reads : int ;
  score : float ;
}
[@@deriving sexp]

val equal : t -> t -> bool

val compare : t -> t -> int

val hash : t -> int

val distance : t -> int option

val contacted_fragment : t -> Genomic_interval.t

val symbol_annotate_baits :
  bait_collection:Genomic_interval_collection.t ->
  genome_annotation:Genomic_annotation.t ->
  max_dist:int ->
  (string list) String.Map.t

val annotations_by_element :
  element_coordinates:Genomic_interval_collection.t ->
  fragments: Genomic_interval_collection.t ->
  fragment_to_baits:(string list) String.Map.t ->
  annotated_baits:(string list) String.Map.t ->
  (string list) String.Map.t

val elements_by_annotation :
  (Genomic_interval.t * string list) list ->
  string list String.Map.t

val output_bait_annotation :
  bait_collection:Genomic_interval_collection.t ->
  bait_annotation:(string list) String.Map.t ->
  path:string ->
  unit

val get_id_frag : t -> string

val get_id_bait : t -> string

val write_contact : t -> Out_channel.t -> unit
