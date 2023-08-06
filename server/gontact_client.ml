open Brr

let () = Console.(log [str "DOM content loaded."])

let decode_mode () =
  Document.body G.document
  |> El.at (Jstr.of_string "data-mode")
  |> Option.map (fun s ->
      Jstr.to_string s
      |> Parsexp.Single.parse_string
      |> Result.map Gontact_shared.mode_of_sexp
    )

let () =
  match decode_mode () with
  | None -> ()
  | Some (Error e) ->
    Console.(log [str (Base.Sexp.to_string (Parsexp.Parse_error.sexp_of_t e))])
  | Some (Ok mode) -> (
      match mode with
      | Form -> ()
      | Results { id = _ } -> ()
    )
