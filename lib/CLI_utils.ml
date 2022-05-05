open Cmdliner


let compare_methods path_annot1 path_annot2 path_output =
    let annot1 = Enhancer_annotation.from_gontact_result_file path_annot1 in
    let annot2 = Enhancer_annotation.from_gontact_result_file path_annot2 in
    Enhancer_annotation.compare_annotations annot1 annot2 path_output

    
let compare_methods_term =
  let open Let_syntax.Cmdliner_term in
  let+ path_annot1 =
    let doc = "Path to GO-enhancer association, first annotation." in
    Arg.(required & opt (some non_dir_file) None & info ["path-annot1"] ~doc ~docv:"PATH")
  and+ path_annot2 =
    let doc = "Path to GO-enhancer association, second annotation." in
    Arg.(required & opt (some non_dir_file) None & info ["path-annot2"] ~doc ~docv:"PATH")
  and+ path_output =
    let doc = "Path to output." in
    Arg.(required & opt (some string) (Some "output.txt") & info ["path-output"] ~doc ~docv:"PATH")
  in
  compare_methods path_annot1 path_annot2 path_output 

let info_compare_methods = Cmd.info ~doc:"GOntact utils compare methods" "compare-methods"
let compare_methods_command = Cmd.v info_compare_methods compare_methods_term

let info = Cmd.info ~doc:"GOntact utils" "gontact-utils"

let command = Cmd.group info [ compare_methods_command ] 
