open Core
open Cmdliner

let main ~mode ~functional_annot ~obo_path ~domain ~gene_info ~fg_path ~bg_path ~bait_coords ~max_dist_bait_TSS ~output_dir ~output_prefix =
  Printf.printf "mode %s functional_annot %s obo_path %s domain %s gene_info %s fg_path %s bg_path %s bait_coords %s max_dist_bait_TSS %d output_dir %s output_prefix %s \n" mode functional_annot obo_path domain gene_info fg_path bg_path bait_coords max_dist_bait_TSS output_dir output_prefix 
  
  (*
  let main_result =
    let open Let_syntax.Result in
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


    *)

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
  and+ bait_coords =
    let doc = "Path to bait coordinates file. " in
    Arg.(required & opt (some non_dir_file) None & info ["bait-coords"] ~doc ~docv:"PATH")
  and+ max_dist_bait_TSS =
    let doc = "Maximum accepted distance (in base pairs) between gene TSS and bait coordinates." in
    Arg.(value & opt int 1_000 & info ["max-dist-bait-TSS"] ~doc ~docv:"INT")
  and+ output_dir =
    let doc = "Output directory." in
    Arg.(required & opt (some string) None & info ["output-dir"] ~doc ~docv:"PATH")
  and+ output_prefix =
    let doc = "Prefix for output files." in
    Arg.(required & opt (some string) None & info ["output-prefix"] ~doc ~docv:"PATH")
   in
  main ~mode ~functional_annot ~obo_path ~domain ~gene_info ~fg_path ~bg_path ~bait_coords ~max_dist_bait_TSS ~output_dir ~output_prefix

let info = Cmd.info ~doc:"Compute GO enrichments." "GOntact"
let command = Cmd.v info term
