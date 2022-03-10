open Core
open Cmdliner

let main ~mode ~functional_annot ~obo_path ~domain ~gene_info ~fg_path ~bg_path ~chr_sizes ~upstream ~downstream ~extend ~bait_coords ~ibed_path ~max_dist_bait_TSS ~output_dir ~output_prefix =
  Printf.printf "mode %s functional_annot %s obo_path %s domain %s gene_info %s fg_path %s bg_path %s bait_coords %s ibed_paths %s max_dist_bait_TSS %d output_dir %s output_prefix %s \n" mode functional_annot obo_path domain gene_info fg_path bg_path bait_coords ibed_path max_dist_bait_TSS output_dir output_prefix  ; 

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
    | "hybrid" -> Error (Printf.sprintf "Hybrid mode not implemented yet.") 
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
        let domains = Great.basal_plus_extension_domains ~chromosome_sizes:chr_collection ~genomic_annotation:filtered_annot_chr ~upstream ~downstream ~extend in
        let domains_int = Great.genomic_interval_collection domains in
        let go_frequencies_foreground = Great.go_frequencies ~element_coordinates:foreground ~regulatory_domains:domains_int ~functional_annot:propagated_fa in
        let go_frequencies_background = Great.go_frequencies ~element_coordinates:background ~regulatory_domains:domains_int ~functional_annot:propagated_fa in
        let enrichment_results = Go_enrichment.foreground_vs_background_binom_test ~go_frequencies_foreground ~go_frequencies_background in
          let output_path = Printf.sprintf "%s/%s_GREAT_results.txt" output_dir output_prefix in
        Go_enrichment.write_output enrichment_results output_path ;
        Ok "GREAT computation finished successfully.";
      )
    | "contacts" -> Error (Printf.sprintf "mode %s not yet implemented.\n" mode) 
    | _ -> Error (Printf.sprintf "mode %s not recognized.\n" mode)
  in
  match main_result with
  | Ok m -> print_endline m
  | Error e -> print_endline e
                 
  
(*
  let output_path = Printf.sprintf "%s/%s.txt"
  let bait_collection = Genomic_interval_collection.of_bed_file bait_coords ~strip_chr:true ~format:Base1 in
      let tss_int = Genomic_annotation.all_tss_intervals filtered_annot_bio_tx max_dist_bait_TSS
    
  in
  let intersection = Genomic_interval_collection.intersect bait_collection tss_int in
  let bait_list = Genomic_interval_collection.interval_list bait_collection in
  let output_one_bait b intersection output =
    let id = Genomic_interval.id b in
    let chr = Genomic_interval.chr b in
    let start_pos = Genomic_interval.start_pos b in
    let end_pos = Genomic_interval.end_pos b in
    let itss = String.Map.find intersection id in
    match itss with
    | None -> Printf.fprintf output "%s\t%d\t%d\t%s\n" chr start_pos end_pos ""
    | Some sl ->
      let ul = List.dedup_and_sort ~compare:String.compare sl in
      Printf.fprintf output "%s\t%d\t%d\t%s\n" chr start_pos end_pos (String.concat ~sep:"," ul)
  in
  Out_channel.with_file output_path ~append:false ~f:(fun output ->
      List.iter bait_list ~f:(fun b -> output_one_bait b intersection output)
    )

*)
    
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
  and+ output_dir =
    let doc = "Output directory." in
    Arg.(value & opt string "." & info ["output-dir"] ~doc ~docv:"PATH")
  and+ output_prefix =
    let doc = "Prefix for output files." in
    Arg.(value & opt string "GOntact" & info ["output-prefix"] ~doc ~docv:"PATH")
   in
  main ~mode ~functional_annot ~obo_path ~domain ~gene_info ~fg_path ~bg_path ~chr_sizes ~upstream ~downstream ~extend ~bait_coords ~ibed_path ~max_dist_bait_TSS ~output_dir ~output_prefix

let info = Cmd.info ~doc:"Compute GO enrichments." "GOntact"
let command = Cmd.v info term
