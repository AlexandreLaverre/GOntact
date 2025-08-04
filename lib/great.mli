open Core

type param = {
  upstream : int ;
  downstream : int ;
  extend : int ;
}

val domains_intervals :
  param ->
  chromosome_sizes:Genomic_interval_collection.t ->
  genome_annotation:Genomic_annotation.t ->
  Genomic_interval_collection.t

type enrichment_analysis = {
  enriched_terms : Go_enrichment.enrichment_result list ;
  domains_int : Genomic_interval_collection.t ;
  element_annotation : Go_enrichment.annotation FGBG.t ;
}

val enrichment_analysis :
  param ->
  chromosome_sizes:Genomic_interval_collection.t ->
  genome_annotation:Genomic_annotation.t ->
  functional_annotation:Functional_annotation.t ->
  Genomic_interval_collection.t FGBG.t ->
  enrichment_analysis

val elements_by_go_category : (Genomic_interval.t * string list) list -> string list String.Map.t

val symbol_elements :
  element_coordinates:Genomic_interval_collection.t ->
  regulatory_domains:Genomic_interval_collection.t ->
  (Genomic_interval.t * Genomic_interval.t list) list

val test_regulatory_domains : unit -> (unit, string) result
