open Core
open Cmdliner

let main ~mode ~functional_annot ~obo_path ~domain ~gene_info ~fg_path ~bg_path ~chr_sizes ~upstream ~downstream ~extend ~bait_coords ~ibed_path ~max_dist_bait_TSS ~max_dist_element_fragment ~min_dist_contacts ~max_dist_contacts ~min_score ~output_dir ~output_prefix ~write_elements_foreground ~write_elements_background =

  let found_func_annot = Sys.file_exists functional_annot in
  let found_obo = Sys.file_exists obo_path in
  let found_gene_info = Sys.file_exists gene_info in
  let found_fg = Sys.file_exists fg_path in
  let found_bg = Sys.file_exists bg_path in
  let found_bait_coords = Sys.file_exists bait_coords in
  let found_chr_sizes = Sys.file_exists chr_sizes in
  let ibed_files = String.split ibed_path ~on:',' in 
  let found_ibed_files = List.for_all ibed_files ~f:(fun file ->
      match (Sys.file_exists file) with
      | `Yes -> true
      | _ -> false
    )
  in
  let check_required_parameters =
    match mode with
    | "GREAT" -> ( 
        match (found_func_annot, found_obo, found_gene_info, found_fg, found_bg, found_chr_sizes) with
        | (`Yes, `Yes, `Yes, `Yes, `Yes, `Yes) -> Ok "All required files were found."
        | _ -> Error "In GREAT mode the following parameters are required: functional-annot, ontology, gene-annot, chr-sizes, foreground, background."
      )
    | "contacts" -> (
        match (found_func_annot, found_obo, found_gene_info, found_fg, found_bg, found_bait_coords, found_ibed_files) with
        | (`Yes, `Yes, `Yes, `Yes, `Yes, `Yes, true) -> Ok "All required files were found."
        | _ -> Error "In contacts mode the following parameters are required: functional-annot, ontology, gene-annot, foreground, background, bait-coords, ibed-path."
      )
    | "hybrid" -> (
        match (found_func_annot, found_obo, found_gene_info, found_fg, found_bg, found_bait_coords, found_chr_sizes, found_ibed_files) with
        | (`Yes, `Yes, `Yes, `Yes, `Yes, `Yes, `Yes, true) -> Ok "All required files were found."
        | _ -> Error "In hybrid mode the following parameters are required: functional-annot, ontology, gene-annot, foreground, background, bait-coords, chr-sizes, ibed-path."
      )
    | x -> Error (Printf.sprintf "Mode %s not recognized." x) 
  in
  let main_result =
    let open Let_syntax.Result in
    let* check_pars = check_required_parameters in print_endline check_pars ; 
    let* namespace = Ontology.define_domain domain in 
    let* obo = Obo.of_obo_file obo_path in
    let* ontology = Ontology.of_obo obo namespace in
    let* gaf = Gaf.of_gaf_file functional_annot in 
    let* gene_annot = Genomic_annotation.of_ensembl_biomart_file gene_info in
    let gonames = Ontology.term_names ontology in 
    let fa = Functional_annotation.of_gaf_and_ontology gaf ontology in
    let propagated_fa = Functional_annotation.propagate_annotations fa ontology in
    let gene_symbols = Functional_annotation.gene_symbols propagated_fa in 
    let filtered_annot_bio_gene = Genomic_annotation.filter_gene_biotypes gene_annot "protein_coding" in (*take only protein_coding genes*)
    let filtered_annot_bio_tx = Genomic_annotation.filter_transcript_biotypes filtered_annot_bio_gene "protein_coding" in (*take only protein_coding transcripts*)
    let filtered_annot = Genomic_annotation.filter_gene_symbols filtered_annot_bio_tx gene_symbols in  (*take only genes whose symbols are in functional (GO) annotations*)
    let foreground = Genomic_interval_collection.of_bed_file fg_path ~strip_chr:true ~format:Base0 in
    let background = Genomic_interval_collection.of_bed_file bg_path ~strip_chr:true ~format:Base0 in
    match mode with
    | "GREAT" -> (
        let chr_collection = Genomic_interval_collection.of_chr_size_file chr_sizes ~strip_chr:true in
        let chr_set = Genomic_interval_collection.chr_set chr_collection in
        let filtered_annot_chr = Genomic_annotation.filter_chromosomes filtered_annot chr_set in (*take only genes on standard chromosomes*)
        let major_isoforms = Genomic_annotation.identify_major_isoforms_symbols filtered_annot_chr in
        let domains = Great.basal_plus_extension_domains ~chromosome_sizes:chr_collection ~genomic_annotation:filtered_annot_chr ~upstream ~downstream ~extend in
        let domains_int = Great.genomic_interval_collection domains in
        let go_elements_foreground = Great.go_elements ~element_coordinates:foreground ~regulatory_domains:domains_int ~functional_annot:propagated_fa in
        let go_elements_background = Great.go_elements ~element_coordinates:background ~regulatory_domains:domains_int ~functional_annot:propagated_fa in
        let go_frequencies_foreground = Go_enrichment.go_frequencies go_elements_foreground in
        let go_frequencies_background = Go_enrichment.go_frequencies go_elements_background in
        let enrichment_results = Go_enrichment.foreground_vs_background_binom_test ~go_frequencies_foreground ~go_frequencies_background in
        let output_path = Printf.sprintf "%s/%s_GREAT_results.txt" output_dir output_prefix in
        Go_enrichment.write_output enrichment_results gonames output_path ;
        if write_elements_foreground then (
          let symbol_elements_foreground = Great.symbol_elements ~element_coordinates:foreground ~regulatory_domains:domains_int in
          let dist_gene_elements_foreground = Genomic_annotation.compute_cis_distances symbol_elements_foreground ~gene_annotation:filtered_annot_chr ~major_isoforms:major_isoforms in 
          let output_path_fg_elements = Printf.sprintf "%s/%s_GREAT_element_gene_association_foreground.txt" output_dir output_prefix in
          Genomic_annotation.write_distance_elements ~dist_elements:dist_gene_elements_foreground output_path_fg_elements ;
        ) ;
        if write_elements_background then (
          let symbol_elements_background = Great.symbol_elements ~element_coordinates:background ~regulatory_domains:domains_int in
          let dist_gene_elements_background = Genomic_annotation.compute_cis_distances symbol_elements_background ~gene_annotation:filtered_annot_chr ~major_isoforms:major_isoforms in 
          let output_path_bg_elements = Printf.sprintf "%s/%s_GREAT_element_gene_association_background.txt" output_dir output_prefix in
          Genomic_annotation.write_distance_elements ~dist_elements:dist_gene_elements_background output_path_bg_elements ;
        ) ;
        Ok "GREAT computation finished successfully.";
      )
    | "contacts" -> (
        let bait_collection = Genomic_interval_collection.of_bed_file bait_coords ~strip_chr:true ~format:Base1 in
        (*let nb_baits = List.length (Genomic_interval_collection.interval_list bait_collection) in
          Printf.printf "Found %d baits.\n" nb_baits ; *) 
        let annotated_baits = Chromatin_contact.go_annotate_baits ~bait_collection ~genome_annotation:filtered_annot ~max_dist:max_dist_bait_TSS ~functional_annot:propagated_fa in
        let output_path_baits = Printf.sprintf "%s/%s_bait_annotation.txt" output_dir output_prefix in
        Chromatin_contact.output_bait_annotation ~bait_collection ~bait_annotation:annotated_baits ~path:output_path_baits ; 
        let contact_list = List.map ibed_files ~f:(fun file -> Chromatin_contact.of_ibed_file file ~strip_chr:true) in
        let with_annotated_baits = List.map contact_list ~f:(fun cc -> Chromatin_contact.remove_unannotated_baits ~contacts:cc ~bait_annotation:annotated_baits) in 
        let cis_contacts = List.map with_annotated_baits ~f:(fun cc -> Chromatin_contact.select_cis cc) in
        let range_contacts = List.map cis_contacts ~f:(fun cc -> Chromatin_contact.select_distance cc ~min_dist:(float_of_int min_dist_contacts) ~max_dist:(float_of_int max_dist_contacts)) in
        let score_contacts = List.map range_contacts ~f:(fun cc -> Chromatin_contact.select_min_score cc ~min_score) in
        let unbaited_contacts = List.map score_contacts ~f:(fun cc -> Chromatin_contact.select_unbaited cc ~bait_collection) in
        let all_contacts = List.dedup_and_sort ~compare:Chromatin_contact.compare (List.join unbaited_contacts) in
        let contacted_fragments = Chromatin_contact.extend_fragments ~contacts:all_contacts ~margin:max_dist_element_fragment in
        (* let nb_contacted_fragments = List.length (Genomic_interval_collection.interval_list contacted_fragments) in
          Printf.printf "Found %d contacted fragments.\n" nb_contacted_fragments;*)
        let fragment_to_baits = Chromatin_contact.fragment_to_baits ~contacts:all_contacts in
        let go_elements_foreground = Chromatin_contact.annot_elements ~element_coordinates:foreground ~fragments:contacted_fragments ~fragment_to_baits ~annotated_baits in
        let go_elements_background = Chromatin_contact.annot_elements ~element_coordinates:background ~fragments:contacted_fragments ~fragment_to_baits ~annotated_baits in
        let go_frequencies_foreground = Go_enrichment.go_frequencies go_elements_foreground in
        let go_frequencies_background = Go_enrichment.go_frequencies go_elements_background in
        let enrichment_results = Go_enrichment.foreground_vs_background_binom_test ~go_frequencies_foreground ~go_frequencies_background in
        let output_path = Printf.sprintf "%s/%s_contacts_results.txt" output_dir output_prefix in
        Go_enrichment.write_output enrichment_results gonames output_path ;
        Ok "GOntact computation finished successfully.";
      )
    | "hybrid" -> (
        let chr_collection = Genomic_interval_collection.of_chr_size_file chr_sizes ~strip_chr:true in
        let chr_set = Genomic_interval_collection.chr_set chr_collection in
        let filtered_annot_chr = Genomic_annotation.filter_chromosomes filtered_annot chr_set in (*take only genes on standard chromosomes*)
        let domains = Great.basal_plus_extension_domains ~chromosome_sizes:chr_collection ~genomic_annotation:filtered_annot_chr ~upstream ~downstream ~extend:0 in
        let domains_int = Great.genomic_interval_collection domains in
        let go_elements_great_foreground = Great.go_elements ~element_coordinates:foreground ~regulatory_domains:domains_int ~functional_annot:propagated_fa in
        let go_elements_great_background = Great.go_elements ~element_coordinates:background ~regulatory_domains:domains_int ~functional_annot:propagated_fa in
        let bait_collection = Genomic_interval_collection.of_bed_file bait_coords ~strip_chr:true ~format:Base1 in
        let annotated_baits = Chromatin_contact.go_annotate_baits ~bait_collection ~genome_annotation:filtered_annot ~max_dist:max_dist_bait_TSS ~functional_annot:propagated_fa in
        let output_path_baits = Printf.sprintf "%s/%s_bait_annotation.txt" output_dir output_prefix in
        Chromatin_contact.output_bait_annotation ~bait_collection ~bait_annotation:annotated_baits ~path:output_path_baits ; 
        let contact_list = List.map ibed_files ~f:(fun file -> Chromatin_contact.of_ibed_file file ~strip_chr:true) in
        let with_annotated_baits = List.map contact_list ~f:(fun cc -> Chromatin_contact.remove_unannotated_baits ~contacts:cc ~bait_annotation:annotated_baits) in 
        let cis_contacts = List.map with_annotated_baits ~f:(fun cc -> Chromatin_contact.select_cis cc) in
        let range_contacts = List.map cis_contacts ~f:(fun cc -> Chromatin_contact.select_distance cc ~min_dist:(float_of_int min_dist_contacts) ~max_dist:(float_of_int max_dist_contacts)) in
        let score_contacts = List.map range_contacts ~f:(fun cc -> Chromatin_contact.select_min_score cc ~min_score) in
        let unbaited_contacts = List.map score_contacts ~f:(fun cc -> Chromatin_contact.select_unbaited cc ~bait_collection) in
        let all_contacts = List.dedup_and_sort ~compare:Chromatin_contact.compare (List.join unbaited_contacts) in
        let contacted_fragments = Chromatin_contact.extend_fragments ~contacts:all_contacts ~margin:max_dist_element_fragment in
        (*let nb_contacted_fragments = List.length (Genomic_interval_collection.interval_list contacted_fragments) in
          Printf.printf "Found %d contacted fragments.\n" nb_contacted_fragments;*)
        let fragment_to_baits = Chromatin_contact.fragment_to_baits ~contacts:all_contacts in
        let go_elements_cc_foreground = Chromatin_contact.annot_elements ~element_coordinates:foreground ~fragments:contacted_fragments ~fragment_to_baits ~annotated_baits in
        let go_elements_cc_background = Chromatin_contact.annot_elements ~element_coordinates:background ~fragments:contacted_fragments ~fragment_to_baits ~annotated_baits in
        let go_elements_foreground = Go_enrichment.combine_go_elements go_elements_great_foreground go_elements_cc_foreground in
        let go_elements_background = Go_enrichment.combine_go_elements go_elements_great_background go_elements_cc_background in 
        let go_frequencies_foreground = Go_enrichment.go_frequencies go_elements_foreground in
        let go_frequencies_background = Go_enrichment.go_frequencies go_elements_background in
        let enrichment_results = Go_enrichment.foreground_vs_background_binom_test ~go_frequencies_foreground ~go_frequencies_background in
        let output_path = Printf.sprintf "%s/%s_hybrid_results.txt" output_dir output_prefix in
        Go_enrichment.write_output enrichment_results gonames output_path ;
        Ok "GOntact hybrid computation finished successfully.";
      )
    | _ -> Error (Printf.sprintf "mode %s not recognized.\n" mode)
  in
  match main_result with
  | Ok m -> print_endline m
  | Error e -> print_endline e
    
    
let term =
  let open Let_syntax.Cmdliner_term in
  let+ mode =
    let doc = "Run mode. Can be 'contacts', 'GREAT' or 'hybrid'" in
    Arg.(value & opt string "contacts" & info ["mode"] ~doc ~docv:"STRING")
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
    let doc = "Gene ontology domain. Can be \"biological_process\", \"molecular_function\" or \"cellular_component\"." in
    Arg.(value & opt string "biological_process" & info ["domain"] ~doc ~docv:"INT")
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
    let doc = "Minimum score for fragments." in
    Arg.(value & opt float 5.0 & info ["min-score"] ~doc ~docv:"FLOAT")
  and+ output_dir =
    let doc = "Output directory." in
    Arg.(value & opt string "." & info ["output-dir"] ~doc ~docv:"PATH")
  and+ output_prefix =
    let doc = "Prefix for output files." in
    Arg.(value & opt string "GOntact" & info ["output-prefix"] ~doc ~docv:"PATH")
  and+ write_elements_foreground =
    let doc = "Write output for gene-element association in foreground." in
    Arg.(value & flag  & info ["write-foreground"] ~doc)
  and+ write_elements_background =
    let doc = "Write output for gene-element association in background." in
    Arg.(value & flag  & info ["write-background"] ~doc)
  in
  main ~mode ~functional_annot ~obo_path ~domain ~gene_info ~fg_path ~bg_path ~chr_sizes ~upstream ~downstream ~extend ~bait_coords ~ibed_path ~max_dist_bait_TSS ~max_dist_element_fragment ~min_dist_contacts ~max_dist_contacts ~min_score ~output_dir ~output_prefix ~write_elements_foreground ~write_elements_background 

let info = Cmd.info ~doc:"Compute GO enrichments." "GOntact"
let command = Cmd.v info term
