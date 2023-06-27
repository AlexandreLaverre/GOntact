open Core

type regulatory_domain

val genomic_interval_collection :
  regulatory_domain list ->
  Genomic_interval_collection.t

val basal_plus_extension_domains_one_chr :
  chr:string ->
  chromosome_size:int ->
  genomic_annotation:Genomic_annotation.t ->
  upstream:int ->
  downstream:int ->
  extend:int ->
  regulatory_domain list

val basal_plus_extension_domains :
  chromosome_sizes:Genomic_interval_collection.t ->
  genomic_annotation:Genomic_annotation.t ->
  upstream:int ->
  downstream:int ->
  extend:int ->
  regulatory_domain list

val go_categories_by_element :
  element_coordinates:Genomic_interval_collection.t ->
  regulatory_domains:Genomic_interval_collection.t ->
  functional_annot:Functional_annotation.t ->
  (Genomic_interval.t * Ontology.PKey.t list) list

val elements_by_go_category : (Genomic_interval.t * string list) list -> string list String.Map.t

val symbol_elements :
  element_coordinates:Genomic_interval_collection.t ->
  regulatory_domains:Genomic_interval_collection.t ->
  (Genomic_interval.t * Genomic_interval.t list) list
