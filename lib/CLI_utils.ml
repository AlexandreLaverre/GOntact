open Core
open Cmdliner


let main mode path_annot1 path_annot2 path_output =
  if String.equal mode "compare_methods" then
    let annot1 = Enhancer_annotation.from_gontact_result_file path_annot1 in
    let annot2 = Enhancer_annotation.from_gontact_result_file path_annot2 in
    Enhancer_annotation.compare_annotations annot1 annot2 path_output
  else ()
    
let term =
  let open Let_syntax.Cmdliner_term in
  let+ mode =
    let doc = "Run mode. Can be 'compare_methods'." in
    Arg.(value & opt string "compare_methods" & info ["mode"] ~doc ~docv:"STRING")
  and+ path_annot1 =
    let doc = "Path to GO-enhancer association, first annotation." in
    Arg.(value & opt string "" & info ["path_annot1"] ~doc ~docv:"PATH")
  and+ path_annot2 =
    let doc = "Path to GO-enhancer association, second annotation." in
    Arg.(value & opt string "" & info ["path_annot2"] ~doc ~docv:"PATH")
  and+ path_output =
    let doc = "Path to output." in
    Arg.(value & opt string "" & info ["path_output"] ~doc ~docv:"PATH")
  in
  main mode path_annot1 path_annot2 path_output 

let info = Cmd.info ~doc:"GOntact utils" "GOntact"
let command = Cmd.v info term
