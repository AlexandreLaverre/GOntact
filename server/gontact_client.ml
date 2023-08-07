open Brr

let decode_mode () =
  Document.body G.document
  |> El.at (Jstr.of_string "data-mode")
  |> Option.map (fun s ->
      Jstr.to_string s
      |> Parsexp.Single.parse_string
      |> Result.map Gontact_shared.mode_of_sexp
    )

module Result_mode = struct

  let decode_enriched_term x =
    let open Base.Option.Monad_infix in
    Jv.find x "go_id" >>| Jv.to_string >>= fun go_id ->
    Jv.find x "go_term" >>| Jv.to_string >>= fun go_term ->
    Jv.find x "enrichment" >>| Jv.to_float >>= fun enrichment ->
    Jv.find x "pval" >>| Jv.to_float >>= fun pval ->
    Jv.find x "fdr" >>| Jv.to_float >>= fun fdr ->
    Some { Gontact_shared.go_id ; go_term ; enrichment ; pval ; fdr }

  let fetch_loop id =
    let open Fut.Result_syntax in
    let open Brr_io.Fetch in

    let decode_response resp =
      let open Fut.Result_syntax in
      let body = Response.as_body resp in
      let+ json = Body.json body in
      match Jv.to_string @@ Jv.Jarray.get json 0 with
      | "In_progress" -> Some Gontact_shared.In_progress
      | "Completed" ->
        let table = Jv.Jarray.get json 1 in
        let n = Jv.Jarray.length table in
        let terms =
          List.init n (fun i -> decode_enriched_term (Jv.Jarray.get table i))
          |> Base.Option.all
        in
        Option.map (fun x -> Gontact_shared.Completed x) terms
      | _ -> None
    in
    let u = Jstr.of_string (Printf.sprintf "/run/%s?format=json" id) in
    let rec loop () =
      let* resp = url u in
      if Response.ok resp then
        let* request_status = decode_response resp in
        match request_status with
        | Some In_progress ->
          Fut.bind (Fut.tick ~ms:10_000) loop
        | Some (Completed enriched_terms) ->
          Fut.return (Ok enriched_terms)
        | None -> Fut.error (Jv.Error.v (Jstr.v "decoding"))
      else
        Fut.return (Error (Jv.Error.v (Response.status_text resp)))
    in
    loop ()

  let set_visibility el ~on:ev =
    let string_of_visibility = function
      | `Visible -> "visibility:visible"
      | `Hidden -> "visibility:hidden"
      | `Collapse -> "visibility:collapse"
    in
    Note_brr.Elr.set_prop
      (El.Prop.jstr (Jstr.v "style"))
      ~on:(Note.E.map (fun x -> Jstr.v (string_of_visibility x)) ev)
      el

  let create_table result_fetch_ev =
    let head = El.thead [
        El.tr [
          El.th [El.txt' "GO term name"] ;
          El.th [El.txt' "Enrichment"] ;
          El.th [El.txt' "p-value"] ;
          El.th [El.txt' "FDR"] ;
        ]
      ] in
    let body = El.tbody [] in
    Note_brr.Elr.set_children body ~on:(
      Note.E.map
        (fun enriched_terms ->
           let row (term : Gontact_shared.enriched_term) =
             El.tr [
               El.td [El.txt' term.go_term] ;
               El.td [El.txt' (Printf.sprintf "%g" term.enrichment)] ;
               El.td [El.txt' (Printf.sprintf "%g" term.pval)] ;
               El.td [El.txt' (Printf.sprintf "%g" term.fdr)] ;
             ]
           in
           List.map row enriched_terms)
        result_fetch_ev
    ) ;
    let at = [At.style (Jstr.v "visibility:collapse")] in
    let t = El.table ~at [head ; body] in
    set_visibility t ~on:(Note.E.map (Fun.const `Visible) result_fetch_ev) ;
    t

  let main id =
    let result_fetch_ev =
      fetch_loop id
      |> Note_brr.Futr.to_event
      |> Note.E.filter_map Base.Result.ok
    in
    let status_bar =
      let p = El.p [El.txt' "Work in progress..."] in
      set_visibility p ~on:(Note.E.map (Fun.const `Collapse) result_fetch_ev) ;
      p
    in
    let download_button =
      let at = [At.style (Jstr.v "visibility:collapse")] in
      let b = El.button ~at [El.txt' "Download full table"] in
      set_visibility b ~on:(Note.E.map (Fun.const `Visible) result_fetch_ev) ;
      b
    in
    Option.iter
      (fun div -> El.set_children div [status_bar ; El.p [ download_button ] ; create_table result_fetch_ev])
      (Document.find_el_by_id G.document (Jstr.v "result-table"))
end

let () =
  match decode_mode () with
  | None -> ()
  | Some (Error e) ->
    Console.(log [str (Base.Sexp.to_string (Parsexp.Parse_error.sexp_of_t e))])
  | Some (Ok mode) -> (
      match mode with
      | Form -> ()
      | Results { id } -> Result_mode.main id
    )
