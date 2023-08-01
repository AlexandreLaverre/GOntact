type param = {
  min_dist : float ;
  max_dist : float ;
  min_score : float ;
}

type annotated_bait_collection = {
  baits : Genomic_interval_collection.t ;
  annotation : Ontology.PKey.t list Core.String.Map.t ;
}

val annotate_baits :
  Genomic_interval_collection.t ->
  genome_annotation:Genomic_annotation.t ->
  functional_annotation:Functional_annotation.t ->
  max_dist_bait_TSS:int ->
  annotated_bait_collection

val aggregate_contact_graphs :
  Chromatin_contact_graph.t list ->
  param ->
  annotated_bait_collection ->
  Chromatin_contact_graph.t

type t = {
  enriched_terms : Go_enrichment.enrichment_result list ;
  element_annotation : Go_enrichment.annotation FGBG.t ;
  contacted_fragments : Genomic_interval_collection.t ;
  fragment_to_baits : string list Core.String.Map.t ;
}

val perform :
  margin:int ->
  annotated_bait_collection ->
  Functional_annotation.t ->
  Chromatin_contact_graph.t ->
  Genomic_interval_collection.t FGBG.t ->
  t
