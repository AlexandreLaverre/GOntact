open Core

(* tests *)

let _test_obo () =
let r = Obo.of_obo_file "data/GeneOntology/go-basic.obo" in
match r with
| Error e -> print_endline e
| Ok tl ->
  let nb = List.length tl in Printf.printf "%d terms in file\n" nb ;
  print_endline (Obo.show tl)


let _test_gaf () =
  let g = Gaf.of_gaf_file "data/GeneOntology/goa_human.gaf" in
  match g with
  | Error e -> print_endline e
  | Ok gl ->
    let nb = List.length gl in Printf.printf "%d gene-GO associations in file\n" nb ;
    print_endline (Gaf.show gl)

let _test_annot_propagation () =
  let open Let_syntax.Result in
  let* obo = Obo.of_obo_file "data/GeneOntology/go-basic.obo" in
  let* bp = Ontology.of_obo obo Biological_process in
  let+ gaf =  Gaf.of_gaf_file "data/GeneOntology/goa_human.gaf" in
  let fa = Functional_annotation.of_gaf_and_ontology gaf bp in
  let fap = Functional_annotation.propagate_annotations fa bp in
  Functional_annotation.write_annotations fa "original_annot.txt" ;
  Functional_annotation.write_annotations fap "propagated_annot.txt" ;
  Printf.printf "%d obo elements \n" (List.length obo) ;
  Printf.printf "%d gaf elements \n" (List.length gaf)

let _test_overlap_coordinates () =
  let baits = Genomic_interval_collection.of_bed_file "data/PCHi-C/human/hg38.baitmap" ~strip_chr:false ~format:Base1 in
  let tss = Genomic_interval_collection.of_bed_file "tss_coords_1000.txt" ~strip_chr:false ~format:Base1 in
  let intersect = Genomic_interval_collection.intersect baits tss in
  List.iter intersect ~f:(fun (bait, tss_list) ->
      Printf.printf "%s\t%s\n" (Genomic_interval.id bait) (String.concat ~sep:"," (List.map tss_list ~f:Genomic_interval.id))
    ) ;
  Genomic_interval_collection.write_output baits "ordered_baits.txt" ~append:false ;
  Genomic_interval_collection.write_output tss "ordered_tss.txt" ~append:false

let _test_merge_coordinates () =
  let cc = Genomic_interval_collection.of_bed_file "test.bed" ~strip_chr:false ~format:Base0 in
  let merged = Genomic_interval_collection.merge_coordinates cc in
  print_endline (Genomic_interval_collection.show merged)

let _test_genomic_annot () =
  let open Let_syntax.Result in
  let+ ga = Genomic_annotation.of_ensembl_biomart_file "data/ensembl_annotations/human/GeneAnnotation_BioMart_Ensembl102_hg38.txt" in
  let filtered1 = Genomic_annotation.filter_gene_biotypes ga "protein_coding" in
  let filtered2 = Genomic_annotation.filter_transcript_biotypes filtered1 "protein_coding" in
  let major_isoforms = Genomic_annotation.identify_major_isoforms filtered2 in
  Map.iteri major_isoforms ~f:(fun ~key:k ~data:d -> Printf.printf "gene %s major isoform %s\n" k d)

let _test_interval_intersection () =
  let cc1 = Genomic_interval_collection.of_bed_file "test1.bed" ~strip_chr:false ~format:Base1 in
  let cc2 = Genomic_interval_collection.of_bed_file "test2.bed" ~strip_chr:false ~format:Base1 in
  let int = Genomic_interval_collection.intersect cc1 cc2 in
  List.iter int ~f:(fun (u, neighboors) ->
      List.iter neighboors ~f:(fun x -> Printf.printf "%s intersects with %s\n" (Genomic_interval.id u) (Genomic_interval.id x))
    )

let _test_fdr () =
  let pvalues = [ ("p1", 0.0015) ; ("p2", 0.3) ; ("p3", 0.1) ; ("p4", 0.05) ; ("p5", 0.005) ; ("p6", 0.001) ; ("p7", 0.005) ] in
  let fdr = Stats.false_discovery_rates pvalues in
  List.iter fdr ~f:(fun (id, pval) -> Printf.printf "id %s fdr %f\n" id pval)

