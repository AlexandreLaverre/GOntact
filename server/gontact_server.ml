open Printf

module H = Tyxml.Html

(* module Form = struct *)
(*   type _ t = *)
(*     | Return : 'a -> 'a t *)
(*     | Map : ('a -> 'b) * 'a t -> 'b t *)
(*     | Both : 'a t * 'b t -> ('a * 'b) t *)
(*     | Checkbox : { name : string ; label : string ; values : string list } -> string t *)

(*   let return x = Return x *)
(*   let map f x = Map (f, x) *)
(*   let checkbox ~name ~label values = Checkbox { label ; name ; values } *)

(*   let rec to_elts : type s. s t -> _ Tyxml.Html.elt list = function *)
(*     | Return _ -> [] *)
(*     | Map (_, x) -> to_elts x *)
(*     | Both (a, b) -> to_elts a @ to_elts b *)
(*     | Checkbox cbx -> *)
(*       Tyxml.Html.[ *)
(*         fieldset *)
(*           ~legend:(legend [txt cbx.label]) *)
(*           (List.map (fun v -> input ~a:[a_name cbx.name ; a_input_type `Radio ; a_value v] ()) cbx.values) ; *)
(*       ] *)
(* end *)

(* module Form_decoding = struct *)
(*   let ( let* ) = Result.bind *)
(*   let ( let+ ) x f = Result.map f x *)
(*   let string args name = *)
(*     match List.assoc_opt name args with *)
(*     | None -> Error (Printf.sprintf "missing arg %s" name) *)
(*     | Some v -> Ok v *)

(*   let run = function *)
(*     | Ok x -> x *)
(*     | Error msg -> Dream.html ~status:`Bad_Request msg *)
(* end *)


module Form = struct
  open Html_types

  type _ t =
    (* | Return : 'a -> ('a, _) t *)
    | Map : ('a -> 'b) * 'a t -> 'b t
    | Both : 'a t * 'b t -> ('a * 'b) t
    | Radio : { name : string ; label : string ;
                encode : 'a -> string ;
                decode : string -> 'a ;
                values : (string * 'a) list } -> 'a t
    | Textarea : { name : string ; label : string } -> string t
    | File : { name : string ; label : string } -> string t

  (* let return x = Return x *)
  (* let map f x = Map (f, x) *)
  let radio ~name ~label ~encode ~decode values = Radio { label ; name ; values ; encode ; decode }
  let textarea ~name ~label = Textarea { name ; label }
  let file ~name ~label = File {name ; label}

  let ( let+ ) x f = Map (f, x)
  let ( and+ ) x y = Both (x, y)

  let rec to_elts : type u. u t -> [> div] Tyxml.Html.elt list = function
    (* | Return _ -> [] *)
    | Map (_, x) -> to_elts x
    | Both (a, b) -> to_elts a @ to_elts b
    | Radio cbx ->
      let open Tyxml.Html in
      let mk_input i (str_repr, v) =
        let id = sprintf "%s-%d" cbx.name i in
        let lab = label ~a:[a_label_for id ; a_class ["form-check-label"]] [txt str_repr] in
        div ~a:[a_class ["form-check"]] [
          input ~a:[a_name cbx.name ; a_id id ; a_input_type `Radio ;
                    a_class ["form-check-input"]; a_value (cbx.encode v)] () ;
          lab ;
        ]
      in
      [ div ~a:[a_class ["mb-3"]] ((label [txt cbx.label]) :: (List.mapi mk_input cbx.values)) ]
    | Textarea { label = l ; name } ->
      let open Tyxml.Html in
      [
        div ~a:[a_class ["mb-3"]] [
          label ~a:[a_label_for name ; a_class ["form-label"]] [txt l] ;
          br () ;
          textarea ~a:[a_name name] (txt "") ;
        ]
      ]
    | File {name ; label = l} ->
      let open Tyxml.Html in
      [
        div ~a:[a_class ["mb-3"]] [
          label ~a:[a_label_for name ; a_class ["form-label"]] [txt l] ;
          input ~a:[a_name name ; a_input_type `File ; a_class ["form-control"]] ()
        ]
      ]

  let to_elt request x =
    let open Tyxml.Html in
    form ~a:[a_method `Post ; a_class ["w-50"]] [
      Unsafe.data @@ Dream.csrf_tag request ;
      div (to_elts x) ;
      input ~a:[a_class ["btn";"btn-primary"]; a_value "Submit"] () ;
    ]
end

let css_link ?a href = H.link ?a ~rel:[`Stylesheet] ~href ()
let ext_js_script ?(a = []) href = H.script ~a:(H.a_src href :: a) (H.txt "")

let html_page ~title body =
  H.html
    (H.head (H.title (H.txt title)) [
        H.meta ~a:[H.a_charset "utf-8"] () ;
        css_link "https://cdn.jsdelivr.net/npm/bootstrap@5.3.0-alpha1/dist/css/bootstrap.min.css" ~a:[H.a_integrity "sha384-GLhlTQ8iRABdZLl6O3oVMWSktQOp6b7In1Zl3/Jr59b6EGGoI1aFkw7cmDA6j6gD" ; H.a_crossorigin `Anonymous] ;
      ])
    (H.body [H.div ~a:[H.a_class ["container"]] body])

module Form_page = struct
  let main_form =
    let open Form in
    let+ species_assembly =
      radio ~name:"species-assembly" ~label:"Species assembly" ~encode:Fun.id ~decode:Fun.id [
        "hg38", "hg38" ;
        "mm10", "mm10" ;
      ]
    and+ test_region_data = textarea ~name:"test-region-data" ~label:"Enter test regions (BED format)"
    and+ test_region_file = file ~name:"test-region-file" ~label:"Or upload a BED file"
    in
    (species_assembly, test_region_data, test_region_file)

  let render request =
    let open Tyxml.Html in
    html_page ~title:"GOntact" [
      h1 [txt "GOntact"] ;
      Form.to_elt request main_form ;
      ext_js_script "https://cdn.jsdelivr.net/npm/bootstrap@5.3.0-alpha1/dist/js/bootstrap.bundle.min.js" ~a:[H.a_integrity "sha384-w76AqPfDkMBDXo30jS1Sgez6pr3x5MlQ1ZAGC+nuZB+EYdgRZgiwxhTBTkF7CXvN" ; H.a_crossorigin `Anonymous] ;
    ]
end

(* let species_assembly_name = "species-assembly" *)

(* let form_page request = *)
(*   let open Tyxml.Html in *)
(*   html_page ~title:"GOntact" [ *)
(*     h1 [txt "GOntact"] ; *)
(*     form ~a:[a_method `Post] [ *)
(*       Unsafe.data @@ Dream.csrf_tag request ; *)
(*       (\* fieldset *\) *)
(*       (\*   ~legend:(legend [txt "Species assembly"]) *\) *)
(*       (\*   [ *\) *)
(*       div ~a:[a_class ["form-check"]] [ *)
(*         input ~a:[a_id "assembly-hg38" ; a_name species_assembly_name ; a_class ["form-check-input"]; a_input_type `Radio ; a_value "hg38" ] () ; *)
(*         label ~a:[a_label_for "assembly-hg38" ; a_class ["form-label"]] [ txt "Human genome" ] ; *)
(*       ] ; *)
(*       div ~a:[a_class ["form-check"]] [ *)
(*         input ~a:[a_id "assembly-mm10" ; a_name species_assembly_name ; a_class ["form-check-input"]; a_input_type `Radio ; a_value "mm10" ] () ; *)
(*         label ~a:[a_label_for "assembly-mm10" ; a_class ["form-label"]] [ txt "Human genome" ] ; *)
(*       ] ; *)
(*       input ~a:[a_input_type `Submit] () ; *)
(*     ] ; *)
(*     ext_js_script "https://cdn.jsdelivr.net/npm/bootstrap@5.3.0-alpha1/dist/js/bootstrap.bundle.min.js" ~a:[H.a_integrity "sha384-w76AqPfDkMBDXo30jS1Sgez6pr3x5MlQ1ZAGC+nuZB+EYdgRZgiwxhTBTkF7CXvN" ; H.a_crossorigin `Anonymous] ; *)
(*   ] *)

let result_page ~species_assembly _request =
  let open Tyxml.Html in
  html_page ~title:"GOntact" [
    txt species_assembly ;
  ]

let html_to_string html =
  Format.asprintf "%a" (Tyxml.Html.pp ()) html

let () =
  Dream.run
  @@ Dream.logger
  @@ Dream.memory_sessions
  @@ Dream.router [

    Dream.get  "/"
      (fun request ->
        Dream.html (html_to_string @@ Form_page.render request));

    Dream.post "/" (fun request ->
        match%lwt Dream.form request with
        | `Ok ["message", _msg ; "species-assembly", species_assembly] ->
          Dream.html (html_to_string (result_page ~species_assembly request))
        | _ ->
          Dream.empty `Bad_Request);

  ]
