type t = {
  great : Great.enrichment_analysis ;
  contacts : Contact_enrichment_analysis.t ;
  enriched_terms : Go_enrichment.enrichment_result list ;
}

val perform :
  Great.param ->
  chromosome_sizes:Genomic_interval_collection.t ->
  genomic_annotation:Genomic_annotation.t ->
  functional_annotation:Functional_annotation.t ->
  margin:int ->
  annotated_baits:Contact_enrichment_analysis.annotated_bait_collection ->
  contact_graph:Chromatin_contact_graph.t ->
  Genomic_interval_collection.t FGBG.t ->
  t
