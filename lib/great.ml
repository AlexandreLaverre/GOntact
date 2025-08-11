open Core

type regulatory_domain = {
  gene_id : string ;
  gene_symbol : string ;
  chr : string ;
  tss_pos : int ;
  start_pos : int ;
  end_pos : int ;
}

let genomic_interval_collection rdl =
  let interval_of_domain rd =
    Genomic_interval.make ~id:rd.gene_symbol rd.chr rd.start_pos rd.end_pos Genomic_interval.Unstranded
  in
  let il = List.map rdl ~f:interval_of_domain in
  Genomic_interval_collection.of_interval_list il

let basal_domain_of_tss (tss:Genomic_interval.t) ~(genomic_annotation:Genomic_annotation.t)  ~upstream:u ~downstream:d ~chromosome_size:cs =
  let (id, chr, tss_pos, strand) = (Genomic_interval.id tss, Genomic_interval.chr tss, Genomic_interval.start_pos tss,  Genomic_interval.strand tss) in  (* start_pos is tss *)
  let gene_symbol = Genomic_annotation.gene_symbol_exn genomic_annotation id in
  match strand with
  | Genomic_interval.Forward ->
    let new_start = max 1 (tss_pos - u) in
    let new_end = min (tss_pos + d - 1) cs in
    { gene_id = id; gene_symbol; chr ; tss_pos ; start_pos = new_start ; end_pos = new_end }
  | Genomic_interval.Reverse ->
    let new_start = max 1 (tss_pos - d + 1) in
    let new_end = min (tss_pos + u) cs in
    { gene_id = id; gene_symbol; chr ; tss_pos ; start_pos = new_start ; end_pos = new_end }
  | Genomic_interval.Unstranded -> invalid_arg "this gene is unstranded!"

let extend_one_domain (d:Genomic_interval.t) ~(genomic_annotation:Genomic_annotation.t) ~left_boundary ~right_boundary ~extend ~upstream ~downstream ~chromosome_size  =
  (*d is a TSS domain domain*)
  let tss = Genomic_interval.start_pos d in
  let chr = Genomic_interval.chr d in
  let id = Genomic_interval.id d in
  let gene_symbol = Genomic_annotation.gene_symbol_exn genomic_annotation id in
  let basal_domain = basal_domain_of_tss d ~genomic_annotation ~upstream ~downstream ~chromosome_size in
  let current_start = basal_domain.start_pos in
  let current_end = basal_domain.end_pos  in
  let new_start = min current_start (max (tss - extend) (left_boundary + 1)) in
  let new_end = max current_end (min (tss + extend - 1) (right_boundary - 1)) in
  { gene_id = id; gene_symbol; chr ; tss_pos = tss ; start_pos = new_start ; end_pos = new_end }

let extend_domains ~genomic_annotation ~ordered_tss ~extend ~upstream ~downstream ~chromosome_size =
  let extend_one_domain = extend_one_domain ~genomic_annotation ~extend ~upstream ~downstream ~chromosome_size in
  let basal_domain_of_tss = basal_domain_of_tss ~genomic_annotation ~upstream ~downstream ~chromosome_size in
  let rec rec_extend tss_list ~acc ~previous =
    match tss_list with
    | [] -> assert false
    | [ h ]  ->  (extend_one_domain h ~left_boundary:(previous.end_pos) ~right_boundary:(chromosome_size + 1)) :: acc
    | i1 :: (i2 :: _ as t) ->
      let b1 = basal_domain_of_tss i1 in
      let b2 = basal_domain_of_tss i2 in
      let ext = extend_one_domain i1 ~left_boundary:(previous.end_pos) ~right_boundary:(b2.start_pos) in
      rec_extend t ~acc:(ext :: acc) ~previous:b1
  in
  match ordered_tss with
  | [] -> []
  | [ d ] -> [ extend_one_domain d ~left_boundary:1 ~right_boundary:(chromosome_size+1) ]
  | i1 :: (i2 :: _ as t) ->
    let b1 = basal_domain_of_tss i1 in
    let b2 = basal_domain_of_tss i2 in
    extend_one_domain i1 ~left_boundary:0 ~right_boundary:(b2.start_pos) :: (rec_extend t ~acc:[] ~previous:b1)

let basal_plus_extension_domains_one_chr  ~chr ~chromosome_size ~genome_annotation:ga ~upstream ~downstream ~extend =
  let chr_set = String.Set.singleton chr in
  let filtered_annot_chr = Genomic_annotation.filter_chromosomes ga chr_set in (*take only genes on only one chromosome *)
  let major_isoforms = Genomic_annotation.identify_major_isoforms filtered_annot_chr in      (*canonical isoform for each gene*)
  let major_tss = Genomic_annotation.major_isoform_tss filtered_annot_chr ~major_isoforms in         (*genomic_interval collection TSS coordinates, they are ordered*)
  let domains_list = extend_domains ~genomic_annotation:ga ~ordered_tss:(Genomic_interval_collection.interval_list major_tss) ~extend ~upstream ~downstream ~chromosome_size in
  domains_list

let basal_plus_extension_domains ~(chromosome_sizes:Genomic_interval_collection.t) ~genome_annotation ~upstream ~downstream ~extend =
  let cl = Genomic_interval_collection.interval_list chromosome_sizes in
  List.concat_map cl ~f:(fun i ->
      basal_plus_extension_domains_one_chr
        ~chr:(Genomic_interval.chr i)
        ~chromosome_size:(Genomic_interval.end_pos i)
        ~genome_annotation ~upstream ~downstream ~extend)

let go_categories_by_element ~(element_coordinates:Genomic_interval_collection.t) ~(regulatory_domains:Genomic_interval_collection.t) ~(functional_annot:Functional_annotation.t) =
  (*regulatory domains were constructed for genes that have at least one GO annotation*)
  (*all their IDs should be in the functional annotation*)
  let intersection = Genomic_interval_collection.intersect element_coordinates regulatory_domains in (*dictionary, element ids -> list of gene symbols*)
  List.Assoc.map intersection ~f:(fun neighboors ->
      let cat_lists =
        List.map neighboors ~f:(fun n ->
            let s = Genomic_interval.id n in
            Functional_annotation.extract_terms_exn functional_annot (`Symbol s)
          )
      in
      GO_term_set.of_sorted_lists_unsafe cat_lists
    )

let elements_by_go_category unique_gocat_by_element =
  let t = String.Table.create () in
  List.iter unique_gocat_by_element ~f:(fun (elt, go_cats) ->
      List.iter go_cats ~f:(fun go_cat -> Hashtbl.add_multi t ~key:go_cat ~data:(Genomic_interval.id elt))
    ) ;
  Hashtbl.to_alist t
  |> String.Map.of_alist_exn

let symbol_elements ~(element_coordinates:Genomic_interval_collection.t) ~(regulatory_domains:Genomic_interval_collection.t) =
   let intersection = Genomic_interval_collection.intersect element_coordinates regulatory_domains in (*dictionary, element ids -> list of gene symbols*)
   intersection

type param = {
  upstream : int ;
  downstream : int ;
  extend : int ;
}

let domains_intervals { upstream ; downstream ; extend } ~chromosome_sizes ~genome_annotation =
  let domains =
    basal_plus_extension_domains
      ~chromosome_sizes ~genome_annotation
      ~upstream ~downstream ~extend
  in
  genomic_interval_collection domains

type enrichment_analysis = {
  enriched_terms : Go_enrichment.enrichment_result list ;
  domains_int : Genomic_interval_collection.t ;
  element_annotation : Go_enrichment.annotation FGBG.t ;
}

let enrichment_analysis
    param ~chromosome_sizes
    ~genome_annotation ~functional_annotation
    elements =
  let domains_int = domains_intervals ~chromosome_sizes ~genome_annotation param in
  let annotate element_coordinates =
    go_categories_by_element
      ~element_coordinates ~regulatory_domains:domains_int
      ~functional_annot:functional_annotation in
  let element_annotation = FGBG.map elements ~f:annotate in
  let enriched_terms = Go_enrichment.binom_test element_annotation functional_annotation in
  { enriched_terms ; domains_int ; element_annotation }
