open Core
open Cmdliner

let main ~bait_coords ~gene_info ~max_dist ~output_path ~pc_only =
  (* Printf.printf "%s %s %d %s\n" bait_coords gene_info max_dist output ; *)
  let bait_collection = Genomic_interval_collection.of_bed_file bait_coords ~strip_chr:true ~format:Base1 in
  let gene_annot = Genomic_annotation.of_ensembl_biomart_file gene_info in
  match gene_annot with
  | Error msg -> Printf.fprintf Stdlib.stderr "Genome annotation is not formatted correctly in %s. \n %s\n" gene_info msg; exit 1
  | Ok ga -> (
      let tss_int = (
        match pc_only with
        | true -> (
            let filtered_annot_bio_gene = Genomic_annotation.filter_gene_biotypes ga "protein_coding" in (*take only protein_coding genes*)
            let filtered_annot_bio_tx = Genomic_annotation.filter_transcript_biotypes filtered_annot_bio_gene "protein_coding" in (*take only protein_coding transcripts*)
            Genomic_annotation.all_tss_intervals filtered_annot_bio_tx max_dist)
        | false -> Genomic_annotation.all_tss_intervals ga max_dist
      )
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
    )
    
let term =
  let open Let_syntax.Cmdliner_term in
  let+ bait_coords =
    let doc = "Path to bait coordinates file. " in
    Arg.(required & opt (some non_dir_file) None & info ["bait-coords"] ~doc ~docv:"PATH")
  and+ gene_info =
    let doc = "Path to gene annotation file." in
    Arg.(required & opt (some non_dir_file) None & info ["gene-annot"] ~doc ~docv:"PATH")
  and+ max_dist =
    let doc = "Maximum accepted distance (in base pairs) between gene TSS and bait coordinates." in
    Arg.(value & opt int 1_000 & info ["max-dist"] ~doc ~docv:"INT")
  and+ pc_only =
    let doc = "Use only protein-coding genes and transcripts in the bait annotations." in
    Arg.(value & opt bool true & info ["pc-only"] ~doc ~docv:"BOOL")
  and+ output_path =
    let doc = "Output path." in
    Arg.(required & opt (some string) None & info ["output"] ~doc ~docv:"PATH")
  in
  main ~bait_coords ~gene_info ~max_dist ~output_path ~pc_only

let info = Cmd.info ~doc:"Annotate baits using gene information extracted from Ensembl BioMart." "annotate-baits"
let command = Cmd.v info term
