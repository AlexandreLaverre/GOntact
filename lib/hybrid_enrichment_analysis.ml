type t = {
  great : Great.enrichment_analysis ;
  contacts : Contact_enrichment_analysis.t ;
  enriched_terms : Go_enrichment.enrichment_result list ;
}

let gocat_assignment_combine a1 a2 =
  FGBG.{
    foreground = Go_enrichment.combine_annotations a1.foreground a2.foreground ;
    background = Go_enrichment.combine_annotations a1.background a2.background ;
  }

let perform
    great_param ~chromosome_sizes ~genome_annotation ~functional_annotation
    ~margin ~annotated_baits ~contact_graph elements =
  let great =
    Great.enrichment_analysis
      great_param
      ~chromosome_sizes ~genome_annotation
      ~functional_annotation
      elements
  in
  let contacts =
    Contact_enrichment_analysis.perform
      ~margin
      annotated_baits
      functional_annotation
      contact_graph
      elements
  in
  let element_annotation = gocat_assignment_combine great.element_annotation contacts.element_annotation in
  let enriched_terms = Go_enrichment.binom_test element_annotation functional_annotation in
  { great ; contacts ; enriched_terms }
