open Core

type annotation = (Genomic_interval.t * GO_term_set.t) list

val combine_annotations : annotation -> annotation -> annotation

type enrichment_result = {
  id : string ;
  observed : float ;
  expected : float ;
  count_foreground : int ;
  total_foreground : int ;
  count_background : int ;
  total_background : int ;
  pval : float ;
  fdr : float ;
}

val binom_test :
  annotation FGBG.t ->
  Functional_annotation.t ->
  enrichment_result list

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
  categories_by_element:annotation ->
  Functional_annotation.t ->
  int String.Map.t
