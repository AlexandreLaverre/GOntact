open Core
open Cmdliner

module Sys = Sys_unix

let compare_methods path_annot1 path_annot2 path_output =
    let annot1 = Enhancer_annotation.from_gontact_result_file path_annot1 in
    let annot2 = Enhancer_annotation.from_gontact_result_file path_annot2 in
    Enhancer_annotation.compare_annotations annot1 annot2 path_output

let shared_contacts ibed_paths min_nb_samples path_output =
  let ibed_files = String.split ibed_paths ~on:',' in
  let found_ibed_files = List.for_all ibed_files ~f:(fun file ->
      match (Sys.file_exists file) with
      | `Yes -> true
      | _ -> false
    )
  in
  if found_ibed_files then
    let table l =
      let module H = Hashtbl.Make(Chromatin_contact) in
      let h = H.create () in
      List.iter l ~f:(fun x -> Hashtbl.add_multi h ~key:x ~data:()) ;
      Hashtbl.to_alist h |> List.map ~f:(fun (x, y) -> (x, List.length y))
    in
    let contact_list = List.concat_map ibed_files ~f:(fun file -> Chromatin_contact_graph.of_ibed_file file ~strip_chr:false) in
    let contact_table = table contact_list in
    Out_channel.with_file path_output ~append:false ~f:(fun output -> (
        Printf.fprintf output "bait_chr\tbait_start\tbait_end\tbait_name\totherEnd_chr\totherEnd_start\totherEnd_end\totherEnd_name\tN_reads\tscore\n" ;
        List.iter contact_table ~f:(fun (c, nb) -> (if nb >= min_nb_samples then Chromatin_contact.write_contact c output))))


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


let shared_contacts_term =
  let open Let_syntax.Cmdliner_term in
  let+ ibed_paths =
    let doc = "Path(s) to chromatin contact data in IBED format. Can be several paths separated by commas. " in
    Arg.(required & opt (some string) None & info ["ibed-path"] ~doc ~docv:"PATH")
  and+ min_nb_samples =
    let doc = "Minimum number of samples in which contacts should be observed." in
    Arg.(value & opt int 1 & info ["min-nb-samples"] ~doc ~docv:"INT")
  and+ path_output =
    let doc = "Path to output." in
    Arg.(required & opt (some string) (Some "output.txt") & info ["path-output"] ~doc ~docv:"PATH")
  in
  shared_contacts ibed_paths min_nb_samples path_output

(* commands *)

let info_compare_methods = Cmd.info ~doc:"GOntact utils compare methods" "compare-methods"
let compare_methods_command = Cmd.v info_compare_methods compare_methods_term

let info_shared_contacts = Cmd.info ~doc:"GOntact utils shared contacts" "shared-contacts"
let shared_contacts_command = Cmd.v info_shared_contacts shared_contacts_term


(* joint command *)

let info = Cmd.info ~doc:"GOntact utils" "gontact-utils"

let command = Cmd.group info [ compare_methods_command ; shared_contacts_command ]
