open Core
open Cmdliner

module Sys = Sys_unix

type parlist = {
  mode : mode ;
  functional_annot : string ;
  obo_path : string ;
  domain : Ontology.domain ;
  gene_info : string ;
  fg_path : string ;
  bg_path : string ;
  chr_sizes : string ;
  upstream : int ;
  downstream : int ;
  extend : int ;
  bait_coords : string ;
  ibed_files : string list ;
  max_dist_bait_TSS : int ;
  max_dist_element_fragment : int ;
  min_dist_contacts : int ;
  max_dist_contacts : int ;
  min_score : float ;
  output_dir : string ;
  output_prefix : string ;
  write_elements_foreground : bool ;
  write_elements_background : bool ;
}
and mode = GREAT | Contacts | Hybrid
[@@deriving sexp]

let save_parlist pl path =
  Sexp.save_hum path (sexp_of_parlist pl)

let check_required_parameters  {mode ;functional_annot ; obo_path ; gene_info ; fg_path ; bg_path ; chr_sizes ; bait_coords ; ibed_files ; _} =
  let found_func_annot = Sys.file_exists functional_annot in
  let found_obo = Sys.file_exists obo_path in
  let found_gene_info = Sys.file_exists gene_info in
  let found_fg = Sys.file_exists fg_path in
  let found_bg = Sys.file_exists bg_path in
  let found_bait_coords = Sys.file_exists bait_coords in
  let found_chr_sizes = Sys.file_exists chr_sizes in
  let found_ibed_files = List.for_all ibed_files ~f:(fun file ->
      match Sys.file_exists file with
      | `Yes -> true
      | _ -> false
    )
  in
  match mode with
  | GREAT -> (
      match found_func_annot, found_obo, found_gene_info, found_fg, found_bg, found_chr_sizes with
      | `Yes, `Yes, `Yes, `Yes, `Yes, `Yes -> Ok "All required files were found."
      | _ -> Error "In GREAT mode the following parameters are required: functional-annot, ontology, gene-annot, chr-sizes, foreground, background."
    )
  | Contacts -> (
      match found_func_annot, found_obo, found_gene_info, found_fg, found_bg, found_bait_coords, found_ibed_files, found_chr_sizes with
      | `Yes, `Yes, `Yes, `Yes, `Yes, `Yes, true, `Yes -> Ok "All required files were found."
        | _ -> Error "In contacts mode the following parameters are required: functional-annot, ontology, gene-annot, chr-sizes, foreground, background, bait-coords, ibed-path."
    )
  | Hybrid -> (
      match found_func_annot, found_obo, found_gene_info, found_fg, found_bg, found_bait_coords, found_chr_sizes, found_ibed_files with
      | `Yes, `Yes, `Yes, `Yes, `Yes, `Yes, `Yes, true -> Ok "All required files were found."
      | _ -> Error "In hybrid mode the following parameters are required: functional-annot, ontology, gene-annot, foreground, background, bait-coords, chr-sizes, ibed-path."
    )

let output_file {output_dir ; output_prefix ; _} suff =
  Printf.sprintf "%s/%s_%s" output_dir output_prefix suff

let expand_go_term_sets go_terms_by_elt fa =
  List.Assoc.map go_terms_by_elt ~f:(fun go_term_set ->
      GO_term_set.to_sorted_list go_term_set
      |> Functional_annotation.term_names_of_pkeys fa
    )

type gocat_assignment = {
  foreground : (Genomic_interval.t * GO_term_set.t) list ;
  background : (Genomic_interval.t * GO_term_set.t) list ;
}

let gocat_assignment_combine a1 a2 =
  {
    foreground = Go_enrichment.combine_GO_maps a1.foreground a2.foreground ;
    background = Go_enrichment.combine_GO_maps a1.background a2.background ;
  }

let enrichment { foreground ; background } propagated_fa =
  let go_frequencies_foreground = Go_enrichment.go_frequencies ~categories_by_element:foreground propagated_fa in
  let go_frequencies_background = Go_enrichment.go_frequencies ~categories_by_element:background propagated_fa in
  Go_enrichment.foreground_vs_background_binom_test ~go_frequencies_foreground ~go_frequencies_background

let great_gocat_assignment ~chromosome_sizes ~genomic_annotation ~upstream ~downstream ~extend ~foreground ~background ~propagated_fa =
  let domains = Great.basal_plus_extension_domains ~chromosome_sizes ~genomic_annotation ~upstream ~downstream ~extend in
  let domains_int = Great.genomic_interval_collection domains in
  let assign element_coordinates =
    Great.go_categories_by_element ~element_coordinates ~regulatory_domains:domains_int
                                   ~functional_annot:propagated_fa in
  { foreground = assign foreground ; background = assign background },
  domains_int

let output_domains pl domains_int =
  let output_path_domains = output_file pl "regulatory_domains.txt" in
  Genomic_interval_collection.write_output domains_int output_path_domains ~append:false

let output_enrichment pl enrichment_results gonames =
  let output_path = output_file pl "enrichment_results.txt" in
  Go_enrichment.write_output enrichment_results gonames output_path

let output_major_isoforms pl major_isoforms =
  let output_path_isoforms = output_file pl "major_isoforms.txt" in
  Genomic_annotation.write_major_isoforms major_isoforms output_path_isoforms

let great_gene_distance_by_element elements ~domains ~filtered_annot ~major_isoforms =
  let element_map = Genomic_interval_collection.interval_map elements in
  let symbol_elements =
    Great.symbol_elements ~element_coordinates:elements ~regulatory_domains:domains
    |> List.map ~f:(fun (elt, neighboors) ->
        Genomic_interval.id elt,
        List.map neighboors ~f:Genomic_interval.id)
    |> String.Map.of_alist_exn
  in
  Genomic_annotation.compute_cis_distances
    symbol_elements
    ~element_map
    ~gene_annotation:filtered_annot
    ~major_isoforms:major_isoforms

let great_mode pl ~chromosome_sizes ~gonames ~filtered_annot ~foreground ~background ~propagated_fa =
  let gocat_assignment, domains_int =
    great_gocat_assignment
      ~chromosome_sizes ~genomic_annotation:filtered_annot
      ~upstream:pl.upstream ~downstream:pl.downstream ~extend:pl.extend
      ~foreground ~background ~propagated_fa
  in
  let enrichment_results = enrichment gocat_assignment propagated_fa in
  output_domains pl domains_int ;
  output_enrichment pl enrichment_results gonames ;

  if (pl.write_elements_foreground || pl.write_elements_background) then (
    let major_isoforms = Utils.chrono "extract major isoforms symbols" Genomic_annotation.identify_major_isoforms_symbols filtered_annot in
    output_major_isoforms pl major_isoforms ;
    if pl.write_elements_foreground then (
      let gocat_by_element_foreground = expand_go_term_sets gocat_assignment.foreground propagated_fa in
      let elements_by_gocat_foreground = Utils.chrono "elements by GO category foreground" Great.elements_by_go_category gocat_by_element_foreground in
      let dist_gene_elements_foreground =
        great_gene_distance_by_element foreground ~domains:domains_int ~filtered_annot ~major_isoforms in
      let output_path_fg_elements = output_file pl "element_gene_association_foreground.txt" in
      Genomic_annotation.write_distance_elements ~dist_elements:dist_gene_elements_foreground output_path_fg_elements ;
      let output_path_fg_elements_go = output_file pl "element_GO_association_foreground.txt" in
      Go_enrichment.write_detailed_association elements_by_gocat_foreground output_path_fg_elements_go ;
    ) ;
    if pl.write_elements_background then (
      let gocat_by_element_background = expand_go_term_sets gocat_assignment.background propagated_fa in
      let elements_by_gocat_background = Utils.chrono "elements by GO category background" Great.elements_by_go_category gocat_by_element_background in
      let dist_gene_elements_background =
        great_gene_distance_by_element background ~domains:domains_int ~filtered_annot ~major_isoforms
      in
      let output_path_bg_elements = output_file pl "element_gene_association_background.txt" in
      Genomic_annotation.write_distance_elements ~dist_elements:dist_gene_elements_background output_path_bg_elements ;
      let output_path_bg_elements_go = output_file pl "element_GO_association_background.txt" in
      Go_enrichment.write_detailed_association elements_by_gocat_background output_path_bg_elements_go ;
    ) ;
  )

let contacts_gocat_assignment ~contacted_fragments ~fragment_to_baits ~annotated_baits ~foreground ~background =
  let assign element_coordinates =
    Chromatin_contact.go_categories_by_element
      ~element_coordinates ~fragments:contacted_fragments
      ~fragment_to_baits ~annotated_baits in
  { foreground = assign foreground ; background = assign background }

let contact_gene_distance_by_element
    elements ~filtered_annot ~major_isoforms
    ~contacted_fragments ~fragment_to_baits ~symbol_annotated_baits =
  let symbol_elements =
    Chromatin_contact.annotations_by_element
      ~element_coordinates:elements
      ~fragments:contacted_fragments
      ~fragment_to_baits
      ~annotated_baits:symbol_annotated_baits
  in
  let element_map = Genomic_interval_collection.interval_map elements in
  Genomic_annotation.compute_cis_distances
    symbol_elements
    ~element_map
    ~gene_annotation:filtered_annot ~major_isoforms

let all_contacts_of_beds beds ~bait_annotation ~bait_collection ~min_dist_contacts ~max_dist_contacts ~min_score =
  let contact_list = List.map beds ~f:(fun file -> Chromatin_contact.of_ibed_file file ~strip_chr:true) in
  let with_annotated_baits = List.map contact_list ~f:(fun cc ->
      Chromatin_contact.remove_unannotated_baits ~contacts:cc ~bait_annotation
    ) in
  let cis_contacts = List.map with_annotated_baits ~f:(fun cc -> Chromatin_contact.select_cis cc) in
  let range_contacts = List.map cis_contacts ~f:(fun cc -> Chromatin_contact.select_distance cc ~min_dist:min_dist_contacts ~max_dist:max_dist_contacts) in
  let score_contacts = List.map range_contacts ~f:(fun cc -> Chromatin_contact.select_min_score cc ~min_score) in
  let unbaited_contacts = List.map score_contacts ~f:(fun cc -> Chromatin_contact.select_unbaited cc ~bait_collection) in
  List.dedup_and_sort ~compare:Chromatin_contact.compare (List.join unbaited_contacts)

let contacts_mode pl ~gonames ~filtered_annot ~foreground ~background ~propagated_fa =
  let bait_collection = Genomic_interval_collection.of_bed_file pl.bait_coords ~strip_chr:true ~format:Base1 in
  let annotated_baits = Chromatin_contact.go_annotate_baits ~bait_collection ~genome_annotation:filtered_annot ~max_dist:pl.max_dist_bait_TSS ~functional_annot:propagated_fa in
  let all_contacts =
    all_contacts_of_beds pl.ibed_files
      ~bait_annotation:annotated_baits ~bait_collection
      ~min_dist_contacts:(float_of_int pl.min_dist_contacts)
      ~max_dist_contacts:(float_of_int pl.max_dist_contacts)
      ~min_score:pl.min_score in
  let contacted_fragments = Chromatin_contact.extend_fragments ~contacts:all_contacts ~margin:pl.max_dist_element_fragment in
  let fragment_to_baits = Chromatin_contact.fragment_to_baits ~contacts:all_contacts in
  let gocat_assignment = contacts_gocat_assignment ~contacted_fragments ~fragment_to_baits ~annotated_baits ~foreground ~background in
  let enrichment_results = enrichment gocat_assignment propagated_fa in
  let annotated_baits = String.Map.map annotated_baits ~f:(Functional_annotation.term_names_of_pkeys propagated_fa) in
  let output_path_GO_baits = output_file pl "bait_GO_annotation.txt" in
  Chromatin_contact.output_bait_annotation ~bait_collection ~bait_annotation:annotated_baits ~path:output_path_GO_baits ;
  output_enrichment pl enrichment_results gonames ;

  if (pl.write_elements_foreground || pl.write_elements_background) then (
    let major_isoforms = Utils.chrono "extract major isoforms symbols" Genomic_annotation.identify_major_isoforms_symbols filtered_annot in
    let symbol_annotated_baits = Chromatin_contact.symbol_annotate_baits ~bait_collection ~genome_annotation:filtered_annot ~max_dist:pl.max_dist_bait_TSS in
    let output_path_gene_baits = output_file pl "bait_gene_annotation.txt" in
    Chromatin_contact.output_bait_annotation ~bait_collection ~bait_annotation:symbol_annotated_baits ~path:output_path_gene_baits ;

    let gene_distance_by_element elements =
      contact_gene_distance_by_element
          elements ~filtered_annot ~major_isoforms
          ~contacted_fragments ~fragment_to_baits ~symbol_annotated_baits
    in
    if pl.write_elements_foreground then (
      let dist_gene_elements_foreground = gene_distance_by_element foreground in
      let output_path_fg_elements = output_file pl "element_gene_association_foreground.txt" in
      Genomic_annotation.write_distance_elements ~dist_elements:dist_gene_elements_foreground output_path_fg_elements ;
      let output_path_fg_elements_go = output_file pl "element_GO_association_foreground.txt" in
      let gocat_by_element_foreground = expand_go_term_sets gocat_assignment.foreground propagated_fa in
      let elements_by_gocat_foreground = Chromatin_contact.elements_by_annotation gocat_by_element_foreground in
      Go_enrichment.write_detailed_association elements_by_gocat_foreground output_path_fg_elements_go ;
    ) ;

    if pl.write_elements_background then (
      let dist_gene_elements_background = gene_distance_by_element background in
      let output_path_bg_elements = output_file pl "element_gene_association_background.txt" in
      Genomic_annotation.write_distance_elements ~dist_elements:dist_gene_elements_background output_path_bg_elements ;
      let output_path_bg_elements_go = output_file pl "element_GO_association_background.txt" in
      let gocat_by_element_background = expand_go_term_sets gocat_assignment.background propagated_fa in
      let elements_by_gocat_background = Chromatin_contact.elements_by_annotation gocat_by_element_background in
      Go_enrichment.write_detailed_association elements_by_gocat_background output_path_bg_elements_go ;
    ) ;
  )

let hybrid_gene_distance_by_element
    elements ~domains ~contacted_fragments ~fragment_to_baits
    ~symbol_annotated_baits ~filtered_annot ~major_isoforms =
  let symbol_elements_cc =
    Chromatin_contact.annotations_by_element
      ~element_coordinates:elements ~fragments:contacted_fragments
      ~fragment_to_baits ~annotated_baits:symbol_annotated_baits in
  let symbol_elements_great =
    Great.symbol_elements ~element_coordinates:elements ~regulatory_domains:domains
    |> List.map ~f:(fun (elt, neighboors) -> Genomic_interval.id elt, List.map neighboors ~f:Genomic_interval.id)
    |> String.Map.of_alist_exn
  in
  let symbol_elements_hybrid = Go_enrichment.combine_maps symbol_elements_cc symbol_elements_great in
  let element_map =  Genomic_interval_collection.interval_map elements in
  Genomic_annotation.compute_cis_distances symbol_elements_hybrid ~element_map ~gene_annotation:filtered_annot ~major_isoforms:major_isoforms

let hybrid_mode pl ~chromosome_sizes ~filtered_annot ~foreground ~background ~propagated_fa ~gonames =
  let gocat_assignment_great, domains_int =
    great_gocat_assignment
      ~chromosome_sizes ~genomic_annotation:filtered_annot
      ~upstream:pl.upstream ~downstream:pl.downstream ~extend:0
      ~foreground ~background ~propagated_fa
  in
  let bait_collection = Genomic_interval_collection.of_bed_file pl.bait_coords ~strip_chr:true ~format:Base1 in
  let annotated_baits = Chromatin_contact.go_annotate_baits ~bait_collection ~genome_annotation:filtered_annot ~max_dist:pl.max_dist_bait_TSS ~functional_annot:propagated_fa in
  let all_contacts =
    all_contacts_of_beds pl.ibed_files
      ~bait_annotation:annotated_baits ~bait_collection
      ~min_dist_contacts:(float_of_int pl.min_dist_contacts)
      ~max_dist_contacts:(float_of_int pl.max_dist_contacts)
      ~min_score:pl.min_score in
  let contacted_fragments = Chromatin_contact.extend_fragments ~contacts:all_contacts ~margin:pl.max_dist_element_fragment in
  let fragment_to_baits = Chromatin_contact.fragment_to_baits ~contacts:all_contacts in
  let gocat_assignment_cc = contacts_gocat_assignment ~annotated_baits ~background ~contacted_fragments ~foreground ~fragment_to_baits in
  let gocat_assignment = gocat_assignment_combine gocat_assignment_great gocat_assignment_cc in
  let enrichment_results = enrichment gocat_assignment propagated_fa in
  let annotated_baits = String.Map.map annotated_baits ~f:(Functional_annotation.term_names_of_pkeys propagated_fa) in
  let output_path_baits = output_file pl "bait_annotation.txt" in
  Chromatin_contact.output_bait_annotation ~bait_collection ~bait_annotation:annotated_baits ~path:output_path_baits ;
  output_enrichment pl enrichment_results gonames ;
  if (pl.write_elements_foreground || pl.write_elements_background) then (
    let major_isoforms = Utils.chrono "extract major isoforms symbols" Genomic_annotation.identify_major_isoforms_symbols filtered_annot in
    let symbol_annotated_baits = Chromatin_contact.symbol_annotate_baits ~bait_collection ~genome_annotation:filtered_annot ~max_dist:pl.max_dist_bait_TSS in

    let gene_distance_by_element elements =
      hybrid_gene_distance_by_element elements
        ~domains:domains_int ~contacted_fragments ~filtered_annot
        ~fragment_to_baits ~symbol_annotated_baits ~major_isoforms in
    if pl.write_elements_foreground then (
      let dist_gene_elements_foreground = gene_distance_by_element foreground in
      let output_path_fg_elements = output_file pl "element_gene_association_foreground.txt" in
      Genomic_annotation.write_distance_elements ~dist_elements:dist_gene_elements_foreground output_path_fg_elements ;
      let output_path_fg_elements_go = output_file pl "element_GO_association_foreground.txt" in
      let gocat_by_element_great_foreground = expand_go_term_sets gocat_assignment_great.foreground propagated_fa in
      let gocat_by_element_cc_foreground = expand_go_term_sets gocat_assignment_great.background propagated_fa in
      let elements_by_gocat_great_foreground = Great.elements_by_go_category gocat_by_element_great_foreground in
      let elements_by_gocat_cc_foreground = Chromatin_contact.elements_by_annotation gocat_by_element_cc_foreground in
      let elements_by_gocat_foreground = Go_enrichment.combine_maps elements_by_gocat_great_foreground elements_by_gocat_cc_foreground in
      Go_enrichment.write_detailed_association elements_by_gocat_foreground output_path_fg_elements_go ;
    ) ;

    if pl.write_elements_background then (
      let dist_gene_elements_background = gene_distance_by_element background in
      let output_path_bg_elements = output_file pl "element_gene_association_background.txt" in
      Genomic_annotation.write_distance_elements ~dist_elements:dist_gene_elements_background output_path_bg_elements ;
      let output_path_bg_elements_go = output_file pl "element_GO_association_background.txt" in
      let gocat_by_element_great_background = expand_go_term_sets gocat_assignment_great.background propagated_fa in
      let gocat_by_element_cc_background = expand_go_term_sets gocat_assignment_cc.background propagated_fa in
      let elements_by_gocat_great_background = Great.elements_by_go_category gocat_by_element_great_background in
      let elements_by_gocat_cc_background = Chromatin_contact.elements_by_annotation gocat_by_element_cc_background in
      let elements_by_gocat_background = Go_enrichment.combine_maps elements_by_gocat_great_background elements_by_gocat_cc_background in
      Go_enrichment.write_detailed_association elements_by_gocat_background output_path_bg_elements_go ;
    ) ;
  )

let main ({mode ;functional_annot ; obo_path ; domain ; gene_info ; fg_path ; bg_path ; chr_sizes ; _} as params) =
  let main_result =
    let open Let_syntax.Result in
    let* check_pars = check_required_parameters params in
    Logs.info (fun m -> m "%s" check_pars) ;
    let* obo = Obo.of_obo_file obo_path in
    let* ontology = Ontology.of_obo obo domain in
    let* gaf = Gaf.of_gaf_file functional_annot in
    let+ gene_annot = Utils.chrono "read genome annotations" Genomic_annotation.of_ensembl_biomart_file gene_info in
    let gonames = Ontology.term_names ontology in
    let fa = Functional_annotation.of_gaf_and_ontology gaf ontology in
    let propagated_fa = Utils.chrono "propagate GO annotations" (Functional_annotation.propagate_annotations fa) ontology in
    let gene_symbols = Functional_annotation.gene_symbols propagated_fa in
    let filtered_annot_bio_gene = Utils.chrono "filter gene biotypes" (Genomic_annotation.filter_gene_biotypes gene_annot) "protein_coding" in (*take only protein_coding genes*)
    Logs.info (fun m -> m "%d genes after filtering gene biotypes" (Genomic_annotation.number_of_genes filtered_annot_bio_gene)) ;
    let filtered_annot_bio_tx = Utils.chrono "filter transcript biotypes" (Genomic_annotation.filter_transcript_biotypes filtered_annot_bio_gene) "protein_coding" in (*take only protein_coding transcripts*)
    Logs.info (fun m -> m "%d genes after filtering transcript biotypes" (Genomic_annotation.number_of_genes filtered_annot_bio_tx)) ;
    let filtered_annot_gene_symbols = Utils.chrono "filter gene symbols" (Genomic_annotation.filter_gene_symbols filtered_annot_bio_tx) gene_symbols in  (*take only genes whose symbols are in functional (GO) annotations*)
    Logs.info (fun m -> m "%d genes after filtering gene symbols" (Genomic_annotation.number_of_genes filtered_annot_gene_symbols)) ;
    let chromosome_sizes = Utils.chrono "construct chr collection" (fun () -> Genomic_interval_collection.of_chr_size_file chr_sizes ~strip_chr:true) () in
    let chr_set = Genomic_interval_collection.chr_set chromosome_sizes in
    let filtered_annot_chr = Utils.chrono "filter standard chromosomes" (Genomic_annotation.filter_chromosomes filtered_annot_gene_symbols) chr_set in (*take only genes on standard chromosomes*)
    Logs.info (fun m -> m "%d genes on standard chromosomes" (Genomic_annotation.number_of_genes filtered_annot_chr)) ;
    let filtered_annot = Utils.chrono "remove duplicated gene symbols" Genomic_annotation.remove_duplicated_gene_symbols filtered_annot_chr in  (*remove duplicated gene symbols*)
    Logs.info (fun m -> m "%d genes after removing duplicated gene symbols" (Genomic_annotation.number_of_genes filtered_annot)) ;
    let unfiltered_foreground = Utils.chrono "read foreground elements" (fun () -> Genomic_interval_collection.of_bed_file fg_path ~strip_chr:true ~format:Base0) () in
    let unfiltered_background = Utils.chrono "read background elements" (fun () -> Genomic_interval_collection.of_bed_file bg_path ~strip_chr:true ~format:Base0) () in
    let foreground = Utils.chrono "removing duplicated foreground elements" Genomic_interval_collection.remove_duplicated_identifiers unfiltered_foreground in
    let background = Utils.chrono "removing duplicated background elements" Genomic_interval_collection.remove_duplicated_identifiers unfiltered_background in
    match mode with
    | GREAT -> great_mode params ~background ~chromosome_sizes ~foreground ~filtered_annot ~gonames ~propagated_fa
    | Contacts -> contacts_mode params ~background ~foreground ~filtered_annot ~gonames ~propagated_fa
    | Hybrid -> hybrid_mode params ~background ~foreground ~chromosome_sizes ~filtered_annot ~gonames ~propagated_fa
  in
  match main_result with
  | Ok () -> ()
  | Error e -> print_endline e

let term =
  let open Let_syntax.Cmdliner_term in
  let+ mode =
    let contacts_label = "contacts" in
    let great_label = "GREAT" in
    let hybrid_label = "hybrid" in
    let modes = [contacts_label, Contacts ; great_label, GREAT ; hybrid_label, Hybrid] in
    let doc = sprintf "Run mode. Can be '%s', '%s' or '%s'" contacts_label great_label hybrid_label in
    Arg.(value & opt (enum modes) Contacts & info ["mode"] ~doc ~docv:"STRING")
  and+ fg_path =
    let doc = "Path to foreground elements, in BED format (0-based, open-ended intervals)." in
    Arg.(required & opt (some non_dir_file) None & info ["foreground"] ~doc ~docv:"PATH")
  and+ bg_path =
    let doc = "Path to background elements, in BED format (0-based, open-ended intervals)." in
    Arg.(required & opt (some non_dir_file) None & info ["background"] ~doc ~docv:"PATH")
  and+ functional_annot =
    let doc = "Path to gene functional annotation in GAF format." in
    Arg.(required & opt (some non_dir_file) None & info ["functional-annot"] ~doc ~docv:"PATH")
  and+ obo_path =
    let doc = "Path to gene ontology definition in OBO format." in
    Arg.(required & opt (some non_dir_file) None & info ["ontology"] ~doc ~docv:"PATH")
  and+ domain =
    let biological_process_label = "biological_process" in
    let molecular_function_label = "molecular_function" in
    let cellular_component_label = "cellular_component" in
    let domains = Ontology.[
        biological_process_label, Biological_process ;
        molecular_function_label, Molecular_function ;
        cellular_component_label, Cellular_component
      ] in
    let doc =
      sprintf "Gene ontology domain. Can be '%s', '%s' or '%s'."
        biological_process_label molecular_function_label cellular_component_label in
    Arg.(value & opt (enum domains) Biological_process & info ["domain"] ~doc ~docv:"INT")
  and+ gene_info =
    let doc = "Path to gene annotation file extracted from Ensembl BioMart." in
    Arg.(required & opt (some non_dir_file) None & info ["gene-annot"] ~doc ~docv:"PATH")
  and+ chr_sizes =
    let doc = "Path to chromosome size files (tab-separated, chr size)." in
    Arg.(value & opt string "NA" & info ["chr-sizes"] ~doc ~docv:"PATH")
  and+ upstream =
    let doc = "Size of basal regulatory domain upstream of TSS. Used in GREAT mode." in
    Arg.(value & opt int 5_000 & info ["upstream"] ~doc ~docv:"INT")
  and+ downstream =
    let doc = "Size of regulatory domain extension. Used in GREAT mode." in
    Arg.(value & opt int 1_000 & info ["downstream"] ~doc ~docv:"INT")
  and+ extend =
    let doc = "Size of basal regulatory domain upstream of TSS. Used in GREAT mode." in
    Arg.(value & opt int 1_000_000 & info ["extend"] ~doc ~docv:"INT")
  and+ bait_coords =
    let doc = "Path to bait coordinates file. " in
    Arg.(value & opt string "NA" & info ["bait-coords"] ~doc ~docv:"PATH")
  and+ ibed_path =
    let doc = "Path(s) to chromatin contact data in IBED format. Can be several paths separated by commas. " in
    Arg.(value & opt string "NA" & info ["ibed-path"] ~doc ~docv:"PATH")
  and+ max_dist_bait_TSS =
    let doc = "Maximum accepted distance (in base pairs) between gene TSS and bait coordinates." in
    Arg.(value & opt int 1_000 & info ["max-dist-bait-TSS"] ~doc ~docv:"INT")
  and+ max_dist_element_fragment =
    let doc = "Maximum accepted distance (in base pairs) between regulatory elements and contacted restriction fragments." in
    Arg.(value & opt int 5_000 & info ["max-dist-element-fragment"] ~doc ~docv:"INT")
  and+ min_dist_contacts =
    let doc = "Minimum accepted distance (in base pairs) between baits and contacted fragments." in
    Arg.(value & opt int 25_000 & info ["min-dist-contacts"] ~doc ~docv:"INT")
  and+ max_dist_contacts =
    let doc = "Maximum accepted distance (in base pairs) between baits and contacted fragments." in
    Arg.(value & opt int 1_000_000 & info ["max-dist-contacts"] ~doc ~docv:"INT")
  and+ min_score =
    let doc = "Minimum score for contacts." in
    Arg.(value & opt float 5.0 & info ["min-score"] ~doc ~docv:"FLOAT")
  and+ output_dir =
    let doc = "Output directory." in
    Arg.(value & opt string "." & info ["output-dir"] ~doc ~docv:"PATH")
  and+ output_prefix =
    let doc = "Prefix for output files." in
    Arg.(value & opt string "GOntact" & info ["output-prefix"] ~doc ~docv:"PATH")
  and+ verbosity =
    let doc = "Verbosity level (info, debug)" in
    Arg.(value & opt (enum Logs.["info", Some Info ; "debug", Some Debug]) None & info ["verbosity"] ~doc)
  and+ write_elements_foreground =
    let doc = "Write output for gene-element association in foreground." in
    Arg.(value & flag  & info ["write-foreground"] ~doc)
  and+ write_elements_background =
    let doc = "Write output for gene-element association in background." in
    Arg.(value & flag  & info ["write-background"] ~doc)
  in
  let ibed_files = String.split ibed_path ~on:',' in
  let pl = {mode ;functional_annot ; obo_path ; domain ; gene_info ; fg_path ; bg_path ; chr_sizes ; upstream ; downstream ; extend ; bait_coords ; ibed_files ; max_dist_bait_TSS ; max_dist_element_fragment ; min_dist_contacts ; max_dist_contacts ; min_score ; output_dir ; output_prefix ; write_elements_foreground ; write_elements_background} in
  let output_pars = Printf.sprintf "%s/%s_parameters.txt" output_dir output_prefix in
  Logs.set_reporter (Logs.format_reporter ());
  Logs.set_level verbosity ;
  save_parlist pl output_pars ;
  main pl

let info = Cmd.info ~doc:"Compute GO enrichments." "GOntact"
let command = Cmd.v info term
