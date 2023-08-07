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

  let main id =
    let result_fetch_ev = Note_brr.Futr.to_event (fetch_loop id) in
    let download_button =
      let at = [ At.disabled ] in
      let b = El.button ~at [El.txt' "Download full table"] in
      Note_brr.Elr.set_prop
        (El.Prop.bool (Jstr.v "disabled"))
        ~on:(Note.E.map (Fun.const false) result_fetch_ev)
        b ;
      b
    in
    Option.iter
      (fun div -> El.set_children div [download_button])
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
