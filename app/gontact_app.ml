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


let () = test_gaf ()
