open Core
    
type t
[@@deriving sexp]
  
val of_ibed_file :
  string ->
  strip_chr:bool->
  t list 

val of_ibed_file_filtered :
  string ->
  strip_chr:bool ->
  min_dist:float ->
  max_dist:float ->
  min_score:float ->
  bait_map:(String.Set.t) ->
  t list 

val select_min_score :
  t list ->
  min_score:float ->
  t list

val select_cis :
  t list ->
  t list

val select_distance :
  t list ->
  min_dist:float ->
  max_dist:float ->
  t list

val select_unbaited :
  t list ->
  bait_collection:Genomic_interval_collection.t ->
  t list

val symbol_annotate_baits :
  bait_collection:Genomic_interval_collection.t ->
  genome_annotation:Genomic_annotation.t ->
  max_dist:int ->
  (string list) String.Map.t

val go_annotate_baits :
  bait_collection:Genomic_interval_collection.t ->
  genome_annotation:Genomic_annotation.t ->
  max_dist:int ->
  functional_annot:Functional_annotation.t ->
  (string list) String.Map.t

val compare : t -> t -> int

val hash : t -> int

val fragment_to_baits :
  contacts:(t list) ->
  (string list) String.Map.t
    
val contacted_fragment_collection :
  contacts:(t list) ->
  Genomic_interval_collection.t

val extend_fragments :
  contacts:(t list) ->
  margin:int ->
  Genomic_interval_collection.t

val annotations_by_element :
  element_coordinates:Genomic_interval_collection.t ->
  fragments: Genomic_interval_collection.t ->
  fragment_to_baits:(string list) String.Map.t ->
  annotated_baits:(string list) String.Map.t ->
  (string list) String.Map.t

val elements_by_annotation :
  (string list) String.Map.t ->
  (string list) String.Map.t 

val output_bait_annotation :
  bait_collection:Genomic_interval_collection.t ->
  bait_annotation:(string list) String.Map.t ->
  path:string ->
  unit

val remove_unannotated_baits :
  contacts:(t list) ->
  bait_annotation:(string list) String.Map.t ->
  (t list)

val get_id_frag : t -> string

val get_id_bait : t -> string

val write_contact : t -> Out_channel.t -> unit

val equal : t -> t -> bool
