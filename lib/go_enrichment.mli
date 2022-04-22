open Core

type enrichment_result 

val foreground_vs_background_binom_test :
  go_frequencies_foreground:(int String.Map.t) ->
  go_frequencies_background:(int String.Map.t) ->
  enrichment_result list 

val write_output :
  enrichment_result list ->
  string String.Map.t -> 
  string ->
  unit

val write_detailed_association :
  (string list) String.Map.t -> 
  string ->
  unit

val combine_maps :
  (string list) String.Map.t ->
  (string list) String.Map.t ->
  (string list) String.Map.t

val go_frequencies :
  categories_by_element:((string list) String.Map.t) ->
  elements_by_category:((string list) String.Map.t) ->
  int String.Map.t 

