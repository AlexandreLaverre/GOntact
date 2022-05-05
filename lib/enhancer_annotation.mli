open Core
    
val from_gontact_result_file :
  string ->
  (string list) String.Map.t 

val compare_annotations :
  (string list) String.Map.t ->
  (string list) String.Map.t ->
  string ->
  unit
