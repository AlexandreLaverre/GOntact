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

  let create_table result_fetch_ev ~fdr_threshold =
    let head = El.thead [
        El.tr [
          El.th [El.txt' "GO term name"] ;
          El.th [El.txt' "Enrichment"] ;
          El.th [El.txt' "p-value"] ;
          El.th [El.txt' "FDR"] ;
        ]
      ] in
    let body = El.tbody [] in
    Note_brr.Elr.def_children body (
      Note.S.l2
        (fun fdr_threshold enriched_terms ->
           let maybe_row (term : Gontact_shared.enriched_term) =
             if term.Gontact_shared.fdr < fdr_threshold then
               Some (
                 El.tr [
                   El.td [El.txt' term.go_term] ;
                   El.td [El.txt' (Printf.sprintf "%g" term.enrichment)] ;
                   El.td [El.txt' (Printf.sprintf "%g" term.pval)] ;
                   El.td [El.txt' (Printf.sprintf "%g" term.fdr)] ;
                 ]
               )
             else None
           in
           List.filter_map maybe_row enriched_terms
        )
        fdr_threshold
        (Note.S.hold [] result_fetch_ev)
    ) ;
    El.table [head ; body]

  let result_widget result_fetch_ev =
    let status_bar =
      let p = El.p [El.txt' "Work in progress..."] in
      set_visibility p ~on:(Note.E.map (Fun.const `Collapse) result_fetch_ev) ;
      p
    in
    let download_button =
      let b = El.button [El.txt' "Download full table"] in
      b
    in
    let fdr_filter_field_id = Jstr.v "fdr-filter-field" in
    let filter_form, fdr_threshold =
      let table_cell_at = At.style (Jstr.v "display:table-cell") in
      let fdr_input = El.input ~at:[At.id fdr_filter_field_id ; At.value (Jstr.v "0.05");table_cell_at] () in
      let fdr_threshold =
        Note_brr.Evr.on_el Ev.keyup (Fun.const ()) fdr_input
        |> Note.E.map (fun  _ -> Jstr.trim @@ El.prop El.Prop.value fdr_input)
        |> Note.E.filter_map (fun s ->
            let x = Jstr.to_float s in
            if Float.is_nan x then None else Some x
          )
        |> Note.S.hold 0.05
      in
      El.fieldset ~at:[At.style (Jstr.v "display:table")] [
        El.legend [El.txt' "Filter"] ;
        El.p ~at:[At.style (Jstr.v"display:table-row")] [
          El.label ~at:[At.for' fdr_filter_field_id; table_cell_at] [El.txt' "FDR threshold"] ;
          fdr_input ;
        ] ;
      ],
      fdr_threshold
    in
    let at = [At.style (Jstr.v "visibility:collapse")] in
    let div = El.div ~at [
        El.div ~at:[At.style (Jstr.v "display:flex;justify-content:space-evenly")] [El.p [ download_button] ; filter_form ] ;
        El.br () ;
        create_table result_fetch_ev ~fdr_threshold
      ]
    in
    set_visibility div ~on:(Note.E.map (Fun.const `Visible) result_fetch_ev) ;
    [ status_bar ; div ]

  let main id =
    let result_fetch_ev =
      fetch_loop id
      |> Note_brr.Futr.to_event
      |> Note.E.filter_map Base.Result.ok
    in
    Option.iter
      (fun div -> El.set_children div (result_widget result_fetch_ev))
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
