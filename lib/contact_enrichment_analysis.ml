open Core

type param = {
  min_dist : int ;
  max_dist : int ;
  min_score : float ;
  min_samples : int option ;
}

type annotated_bait_collection = {
  baits : Genomic_interval_collection.t ;
  annotation : Ontology.PKey.t list String.Map.t ;
}

let annotate_baits baits ~genome_annotation ~functional_annotation ~max_dist_bait_TSS =
  let tss_intervals = Genomic_annotation.all_tss_intervals genome_annotation max_dist_bait_TSS in
  let intersection = (*String.Map - key = bait ids ; values = list of gene_id:gene_symbol mixed id*)
    Genomic_interval_collection.intersect baits tss_intervals
    |> List.map ~f:(fun (elt, fragments) -> Genomic_interval.id elt, List.map fragments ~f:Genomic_interval.id)
    |> String.Map.of_alist_exn
  in
  let get_terms_symbol id =
    match String.split id ~on:':' with
    | [ _ ; symbol ] -> (
      let terms =  Functional_annotation.extract_terms functional_annotation (`Symbol symbol) in
      match terms with
      | Some l -> l
      | None -> []
    )
    | _ -> []
  in
  let go_annot = String.Map.map intersection ~f:(fun l -> List.concat_map l ~f:get_terms_symbol) in
  let annotation = String.Map.map go_annot ~f:(fun l -> List.dedup_and_sort ~compare:Ontology.PKey.compare l) in
  { baits ; annotation }

type t = {
  enriched_terms : Go_enrichment.enrichment_result list ;
  element_annotation : Go_enrichment.annotation FGBG.t ;
  contacted_fragments : Genomic_interval_collection.t ;
  fragment_to_baits : string list String.Map.t ;
}

let aggregate_contact_graphs graphs { min_dist ; max_dist ; min_score ; min_samples } annotated_baits =
  let transform_contact_graph graph =
    graph
    |> Chromatin_contact_graph.remove_unannotated_baits ~bait_annotation:annotated_baits.annotation
    |> Chromatin_contact_graph.select_cis
    |> Chromatin_contact_graph.select_distance ~min_dist ~max_dist
    |> Chromatin_contact_graph.select_min_score ~min_score
    |> Chromatin_contact_graph.select_unbaited ~bait_collection:annotated_baits.baits
  in
  let select = match min_samples with
    | None -> List.hd
    | Some min_samples -> (fun xs -> if List.length xs >= min_samples then List.hd xs else None)
  in
  List.concat_map graphs ~f:transform_contact_graph
  |> List.sort_and_group ~compare:Chromatin_contact.compare
  |> List.filter_map ~f:select

let go_categories_by_element
    ~(element_coordinates:Genomic_interval_collection.t)
    ~(fragments:Genomic_interval_collection.t)
    ~fragment_to_baits ~bait_annotation =
  (* for each element, find the list of annotations - gene symbols or GO categories - associated with it *)
  let intersection =
    Genomic_interval_collection.intersect element_coordinates fragments
  in
  let elbaits =
    List.Assoc.map intersection ~f:(fun l ->
        List.filter_map l ~f:(fun frag -> Map.find fragment_to_baits (Genomic_interval.id frag))
        |> List.join
        |> List.dedup_and_sort ~compare:String.compare
      )
  in
  List.Assoc.map elbaits ~f:(fun l ->
      List.filter_map l ~f:(fun bait -> Map.find bait_annotation bait)
      |> GO_term_set.of_sorted_lists_unsafe
    )

let perform ~margin annotated_baits functional_annotation contact_graph elements =
  let contacted_fragments = Chromatin_contact_graph.extend_fragments contact_graph ~margin in
  let fragment_to_baits = Chromatin_contact_graph.fragment_to_baits contact_graph in
  let annotate element_coordinates =
    go_categories_by_element
      ~element_coordinates ~fragments:contacted_fragments
      ~fragment_to_baits ~bait_annotation:annotated_baits.annotation in
  let element_annotation = FGBG.map elements ~f:annotate in
  let enriched_terms = Go_enrichment.binom_test element_annotation functional_annotation in
  { enriched_terms ; element_annotation ; contacted_fragments ; fragment_to_baits }
