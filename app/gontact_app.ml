open Core
open Gontact    

let () = 
let r = Obo.from_file "/home/ubuntu/data/mydatalocal/GOntact/data/GeneOntology/go-basic.obo" in
match r with
| Error e -> print_endline e
| Ok tl ->
  let nb = List.length tl in Printf.printf "%d terms in file\n" nb ;
  print_endline (Obo.show tl)
(* List.iter tl ~f:(fun t -> List.iter t.is_a ~f:(Printf.printf "%s %s %s \n" t.id t.name t.namespace) ;*)
  
