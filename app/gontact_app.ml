open Core
open Gontact    

let test_obo () = 
let r = Obo.of_obo_file "/home/ubuntu/data/mydatalocal/GOntact/data/GeneOntology/go-basic.obo" in
match r with
| Error e -> print_endline e
| Ok tl ->
  let nb = List.length tl in Printf.printf "%d terms in file\n" nb ;
  print_endline (Obo.show tl)
  

let test_gaf () =
  let g = Gaf.of_gaf_file "/home/ubuntu/data/mydatalocal/GOntact/data/GeneOntology/goa_human.gaf" in
  match g with
  | Error e -> print_endline e
  | Ok gl ->
    let nb = List.length gl in Printf.printf "%d gene-GO associations in file\n" nb ;
    print_endline (Gaf.show gl)

(*let test_func_annot () =
  let g = Gaf.from_file "/home/ubuntu/data/mydatalocal/GOntact/data/GeneOntology/goa_human.gaf" in
  match g with
  | Error e -> print_endline e
  | Ok gl ->
    let fa = Functional_annotation.of_gaf gl in
    print_endline (Functional_annotation.show fa `Gene_symbol_to_go) 
*)    

let test_annot_propagation () =
  let open Let_syntax.Result in
  let* obo = Obo.of_obo_file "/home/ubuntu/data/mydatalocal/GOntact/data/GeneOntology/go-basic.obo" in 
  let* bp = Ontology.of_obo obo `Biological_process in 
  let+ gaf =  Gaf.of_gaf_file "/home/ubuntu/data/mydatalocal/GOntact/data/GeneOntology/goa_human.gaf" in
  let fa = Functional_annotation.of_gaf_and_ontology gaf bp in
  let fap = Functional_annotation.propagate_annotations fa bp in 
  (*let gl = Functional_annotation.extract_genes fa ~go_id:"GO:0000151" `Symbol in
    print_endline ([%show: string list option] gl) *) 
  Functional_annotation.write_annotations fa "original_annot.txt" ;
  Functional_annotation.write_annotations fap "propagated_annot.txt" ;
  Printf.printf "%d obo elements \n" (List.length obo) ;
  Printf.printf "%d gaf elements \n" (List.length gaf)


let test_merge_coordinates () =
  let cc = Genomic_interval_collection.of_bed_file "test.bed" ~strip_chr:false in
  let merged = Genomic_interval_collection.merge_coordinates cc in
  print_endline (Genomic_interval_collection.show merged)

let test_genomic_annot () =
  let open Let_syntax.Result in
  let+ ga = Genomic_annotation.of_ensembl_biomart_file "/home/ubuntu/data/mydatalocal/GOntact/data/ensembl_annotations/human/GeneAnnotation_BioMart_Ensembl102_hg38.txt" in
  let filtered1 = Genomic_annotation.filter_gene_biotypes ga "protein_coding" in
  let filtered2 = Genomic_annotation.filter_transcript_biotypes filtered1 "protein_coding" in
  let major_isoforms = Genomic_annotation.identify_major_isoforms filtered2 in
  String.Map.iteri major_isoforms ~f:(fun ~key:k ~data:d -> Printf.printf "gene %s major isoform %s\n" k d)

let test_interval_intersection () =
  let cc1 = Genomic_interval_collection.of_bed_file "test1.bed" ~strip_chr:false in
  let cc2 = Genomic_interval_collection.of_bed_file "test2.bed" ~strip_chr:false in
  let int = Genomic_interval_collection.intersect cc1 cc2 in
  String.Map.iteri int ~f:(fun ~key:k ~data:d -> (List.iter d ~f:(fun x -> Printf.printf "%s intersects with %s\n" k x)))

let test_regulatory_domains () = 
  let open Let_syntax.Result in
  let* ga = Genomic_annotation.of_ensembl_biomart_file "/home/ubuntu/data/mydatalocal/GOntact/data/ensembl_annotations/human/GeneAnnotation_BioMart_Ensembl102_hg38.txt" in
  let* obo = Obo.of_obo_file "/home/ubuntu/data/mydatalocal/GOntact/data/GeneOntology/go-basic.obo" in 
  let* bp = Ontology.of_obo obo `Biological_process in 
  let+ gaf =  Gaf.of_gaf_file "/home/ubuntu/data/mydatalocal/GOntact/data/GeneOntology/goa_human.gaf" in
  let fa = Functional_annotation.of_gaf_and_ontology gaf bp in
  let fap = Functional_annotation.propagate_annotations fa bp in
  let gene_symbols = Functional_annotation.gene_symbols fap in
  let chr_sizes = Genomic_interval_collection.of_chr_size_file "/home/ubuntu/data/mydatalocal/GOntact/data/ensembl_annotations/human/chr_sizes_hg38.txt" ~strip_chr:false in
  let chr_set = Genomic_interval_collection.chr_set chr_sizes in
  let filtered_annot_chr = Genomic_annotation.filter_chromosomes ga chr_set in (*take only genes on standard chromosomes*)
  let filtered_annot_bio_gene = Genomic_annotation.filter_gene_biotypes filtered_annot_chr "protein_coding" in (*take only protein_coding genes*)
  let filtered_annot_bio_tx = Genomic_annotation.filter_transcript_biotypes filtered_annot_bio_gene "protein_coding" in (*take only protein_coding transcripts*)
  let filtered_annot = Genomic_annotation.filter_gene_symbols filtered_annot_bio_tx gene_symbols in  (*take only genes whose symbols are in functional (GO) annotations*)
  let domains = Great.basal_plus_extension_domains ~chromosome_sizes:chr_sizes ~genomic_annotation:filtered_annot ~upstream:5000 ~downstream:1000 ~extend:1000000 in
  let int_domains = Great.genomic_interval_collection domains in
  Genomic_interval_collection.write_output int_domains "basal_plus_extension_domains.txt" ~append:false 


let test_go_frequencies () =
  let open Let_syntax.Result in
  let* ga = Genomic_annotation.of_ensembl_biomart_file "/home/ubuntu/data/mydatalocal/GOntact/data/ensembl_annotations/human/GeneAnnotation_BioMart_Ensembl102_hg38.txt" in
  let* obo = Obo.of_obo_file "/home/ubuntu/data/mydatalocal/GOntact/data/GeneOntology/go-basic.obo" in 
  let* bp = Ontology.of_obo obo `Biological_process in 
  let+ gaf =  Gaf.of_gaf_file "/home/ubuntu/data/mydatalocal/GOntact/data/GeneOntology/goa_human.gaf" in
  let fa = Functional_annotation.of_gaf_and_ontology gaf bp in
  let fap = Functional_annotation.propagate_annotations fa bp in
  let gene_symbols = Functional_annotation.gene_symbols fap in
  let chr_sizes = Genomic_interval_collection.of_chr_size_file "/home/ubuntu/data/mydatalocal/GOntact/data/ensembl_annotations/human/chr_sizes_hg38.txt" ~strip_chr:false in
  let chr_set = Genomic_interval_collection.chr_set chr_sizes in
  let filtered_annot_chr = Genomic_annotation.filter_chromosomes ga chr_set in (*take only genes on standard chromosomes*)
  let filtered_annot_bio_gene = Genomic_annotation.filter_gene_biotypes filtered_annot_chr "protein_coding" in (*take only protein_coding genes*)
  let filtered_annot_bio_tx = Genomic_annotation.filter_transcript_biotypes filtered_annot_bio_gene "protein_coding" in (*take only protein_coding transcripts*)
  let filtered_annot = Genomic_annotation.filter_gene_symbols filtered_annot_bio_tx gene_symbols in  (*take only genes whose symbols are in functional (GO) annotations*)
  let domains = Great.basal_plus_extension_domains ~chromosome_sizes:chr_sizes ~genomic_annotation:filtered_annot ~upstream:5000 ~downstream:1000 ~extend:1000000 in
  let domains_int = Great.genomic_interval_collection domains in 
  let fg = Genomic_interval_collection.of_bed_file "/home/ubuntu/data/mydatalocal/GOntact/data/enhancers/human/FANTOM5.first1000.kidney.enhancers.hg38.bed" ~strip_chr:true in
  let bg = Genomic_interval_collection.of_bed_file "/home/ubuntu/data/mydatalocal/GOntact/data/enhancers/human/FANTOM5.Laverre2022.bed" ~strip_chr:true in
  let go_freq_fg = Great.go_frequencies ~element_coordinates:fg ~regulatory_domains:domains_int ~functional_annot:fap in
  let go_freq_bg = Great.go_frequencies ~element_coordinates:bg ~regulatory_domains:domains_int ~functional_annot:fap in
  Great.write_go_frequencies ~go_frequencies:go_freq_fg "GOFrequencies_BiologicalProcess_Foreground.txt" ;
  Great.write_go_frequencies ~go_frequencies:go_freq_bg "GOFrequencies_BiologicalProcess_Background.txt" 

let res = test_go_frequencies () ;
  
