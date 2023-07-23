type t = Chromatin_contact.t list

val of_ibed_file : string -> strip_chr:bool-> t

val of_ibed_file_filtered :
  string ->
  strip_chr:bool ->
  min_dist:float ->
  max_dist:float ->
  min_score:float ->
  bait_map:Core.String.Set.t ->
  t

val select_min_score :
  t ->
  min_score:float ->
  t

val select_cis : t -> t

val select_distance :
  t ->
  min_dist:float ->
  max_dist:float ->
  t

val select_unbaited :
  t ->
  bait_collection:Genomic_interval_collection.t ->
  t

val remove_unannotated_baits :
  t ->
  bait_annotation:_ Core.String.Map.t ->
  t

val extend_fragments :
  t ->
  margin:int ->
  Genomic_interval_collection.t

val fragment_to_baits : t ->  string list Core.String.Map.t

val contacted_fragment_collection : t -> Genomic_interval_collection.t
