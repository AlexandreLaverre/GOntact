open Cmdliner

let main ~bait_coords ~gene_info ~margin =
  Printf.printf "%s %s %d\n" bait_coords gene_info margin

let term =
  let open Let_syntax.Cmdliner_term in
  let+ bait_coords =
    let doc = "Path to bait coords file" in
    Arg.(required & opt (some non_dir_file) None & info ["bait-coords"] ~doc ~docv:"PATH")
  and+ gene_info =
    let doc = "Path to gene annotation" in
    Arg.(required & opt (some non_dir_file) None & info ["gene-annot"] ~doc ~docv:"PATH")
  and+ margin =
    let doc = "Margin (in bp)" in
    Arg.(value & opt int 1_000 & info ["margin"] ~doc ~docv:"INT")
  in
  main ~bait_coords ~gene_info ~margin

let info = Cmd.info ~doc:"Annotate baits" "annotate-baits"
let command = Cmd.v info term
