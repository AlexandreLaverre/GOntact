open Core
open Gontact    

let test_obo () = 
let r = Obo.from_file "/home/ubuntu/data/mydatalocal/GOntact/data/GeneOntology/go-basic.obo" in
match r with
| Error e -> print_endline e
| Ok tl ->
  let nb = List.length tl in Printf.printf "%d terms in file\n" nb ;
  print_endline (Obo.show tl)
  

let test_gaf () =
  let g = Gaf.from_file "/home/ubuntu/data/mydatalocal/GOntact/data/GeneOntology/goa_human.gaf" in
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
    let fa = Functional_annotation.from_gaf gl in
    print_endline (Functional_annotation.show fa `Gene_symbol_to_go) 
*)    

let test_annot_propagation () =
  let open Let_syntax.Result in
  let* obo = Obo.from_file "/home/ubuntu/data/mydatalocal/GOntact/data/GeneOntology/go-basic.obo" in 
  let* bp = Ontology.from_obo obo `Biological_process in 
  let+ gaf =  Gaf.from_file "/home/ubuntu/data/mydatalocal/GOntact/data/GeneOntology/goa_human.gaf" in
  let fa = Functional_annotation.from_gaf_and_ontology gaf bp in
  let fap = Functional_annotation.propagate_annotations fa bp in 
  (*let gl = Functional_annotation.extract_genes fa ~go_id:"GO:0000151" `Symbol in
    print_endline ([%show: string list option] gl) *) 
  Functional_annotation.write_annotations fa "original_annot.txt" ;
  Functional_annotation.write_annotations fap "propagated_annot.txt" ;
  Printf.printf "%d obo elements \n" (List.length obo) ;
  Printf.printf "%d gaf elements \n" (List.length gaf)


let test_merge_coordinates () =
  let cc = Genomic_interval_collection.from_bed_file "test.bed" in
  let merged = Genomic_interval_collection.merge_coordinates cc in
  print_endline (Genomic_interval_collection.show merged)


let res  = test_merge_coordinates () ;
