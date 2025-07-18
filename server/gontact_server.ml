open Tyxml
open Core
open Gontact
open Gontact_shared
open Ppx_yojson_conv_lib.Yojson_conv.Primitives

let css_link ?a href = Html.link ?a ~rel:[`Stylesheet] ~href ()
(* let ext_js_script ?(a = []) href = Html.script ~a:(Html.a_src href :: a) (Html.txt "") *)

let great_route = "/great"

let html_page ?mode ~title:page_title contents =
  let open Html in
  let page_head =
    head (title (txt page_title)) [
      meta ~a:[a_charset "utf-8"] () ;
      meta ~a:[a_http_equiv "X-UA-Compatible" ; a_content "IE=edge"] () ;
      meta ~a:[a_name "viewport" ; a_content "width=device-width, initial-scale=1"] () ;
      link ~rel:[`Icon] ~href:"/img/favicon.ico" () ;
      css_link "/css/marx.css" ;
      css_link "/css/gontact.css" ;
      script ~a:[a_src "/js/gontact_client.js" ; a_defer () ] (txt "") ;
    ]
  in
  let atts = match mode with
    | None -> []
    | Some m -> [a_user_data "mode" (Sexp.to_string_mach (sexp_of_mode m))]
  in
  html
    ~a:[a_lang "en"]
    page_head
    (body ~a:atts [main contents])

let logo_header =
  let open Html in
  header [
    img ~src:"/img/GOntact_logo.png" ~alt:"GOntact logo" ~a:[a_width 600] () ;
  ]

let ibed_directory = function
  | `human -> "data/PCHi-C/human/ibed_files"
  | `mouse -> "data/PCHi-C/mouse/ibed_files"

let ibed_files genome =
  let dir = ibed_directory genome in
  Sys_unix.readdir dir
  |> Array.filter ~f:(Core.String.is_suffix ~suffix:".ibed")
  |> Array.map ~f:(Filename.concat dir)
  |> Array.to_list

let load_genome_annotation genome ~gene_symbols =
  let path = match genome with
    | `human -> "data/ensembl_annotations/human/GeneAnnotation_BioMart_Ensembl102_hg38.txt"
    | `mouse -> "data/ensembl_annotations/mouse/GeneAnnotation_BioMart_Ensembl102_mm10.txt"
  in
  let gene_annot =
    Genomic_annotation.of_ensembl_biomart_file path
    |> Result.ok_or_failwith
  in
  let filtered_annot_bio_gene = Genomic_annotation.filter_gene_biotypes gene_annot "protein_coding" in
  let filtered_annot_bio_tx = Genomic_annotation.filter_transcript_biotypes filtered_annot_bio_gene "protein_coding" in
  Genomic_annotation.filter_gene_symbols filtered_annot_bio_tx gene_symbols

let load_annotated_baits genome ~genome_annotation ~functional_annotation =
  let bait_path = match genome with
    | `human -> "data/PCHi-C/human/hg38.baitmap"
    | `mouse -> "data/PCHi-C/mouse/mm10.baitmap"
  in
  let bait_collection =
    Genomic_interval_collection.of_bed_file bait_path ~strip_chr:true ~format:Base1
  in
  Contact_enrichment_analysis.annotate_baits
    bait_collection
    ~genome_annotation ~functional_annotation
    ~max_dist_bait_TSS:1_000

let load_contact_graph genome ~cea_param ~annotated_baits =
  let contact_graphs =
    ibed_files genome
    |> List.map ~f:(Gontact.Chromatin_contact_graph.of_ibed_file ~strip_chr:true)
  in
  Contact_enrichment_analysis.aggregate_contact_graphs
    contact_graphs
    cea_param
    annotated_baits

let load_chromosome_sizes genome =
  let path = match genome with
    | `human -> "data/ensembl_annotations/human/chr_sizes_hg38.txt"
    | `mouse -> "data/ensembl_annotations/mouse/chr_sizes_mm10.txt"
  in
  Genomic_interval_collection.of_chr_size_file path ~strip_chr:true

let load_functional_annotation genome ~domain =
  (
    let open Let_syntax.Result in
    let* obo = Obo.of_obo_file "data/GeneOntology/go-basic.obo" in
    let* ontology = Ontology.of_obo obo domain in
    let gaf_path = match genome with
      | `human -> "data/GeneOntology/goa_human.gaf"
      | `mouse -> "data/GeneOntology/mgi.gaf"
    in
    let+ gaf = Gaf.of_gaf_file gaf_path in
    let fa = Functional_annotation.of_gaf_and_ontology gaf ontology in
    let propagated_fa = Functional_annotation.propagate_annotations fa ontology in
    propagated_fa, ontology
  )
  |> Result.ok_or_failwith

let load_bed bed_contents =
  String.split_lines bed_contents
  |> Genomic_interval_collection.of_bed_lines ~strip_chr:true ~format:Genomic_interval_collection.Base0

let load_default_background genome =
  In_channel.read_all (
    match genome with
    | `human -> "data/enhancers/human/ENCODE.Laverre2022.bed"
    | `mouse -> "data/enhancers/mouse/ENCODE.Laverre2022.bed"
  )

let request_table = String.Table.create ()

module Decode_form = struct
  let scalar ~of_string_exn label = function
    | [None, v] -> (
        try Ok (of_string_exn v)
        with Failure _ -> Error (`Param_parsing label)
      )
    | _ -> Rresult.R.error_msgf "wrong field structure (%s)" label

  let float = scalar ~of_string_exn:Core.Float.of_string
  let int = scalar ~of_string_exn:Core.Int.of_string

  let genome = scalar ~of_string_exn:(function
      | "human" -> `human
      | "mouse" -> `mouse
      | _ -> failwith "unknown genome"
    )

  let domain = scalar ~of_string_exn:(function
      | "cellular_component" -> Ontology.Cellular_component
      | "molecular_function" -> Molecular_function
      | "biological_process" -> Biological_process
      | _ -> failwith "unknown GO domain"
    )

  let single_file_contents label = function
    | [name, contents] -> Ok (name, contents)
    | _ -> Rresult.R.error_msgf "wrong field structure (%s)" label

  let maybe_single_file_contents label = function
    | [] -> Ok None
    | [name, contents] -> Ok (Some (name, contents))
    | _ -> Rresult.R.error_msgf "wrong field structure (%s)" label

end

module Form_page = struct
  let%html intro_par = {|
    <p style="text-align:justify;">
      GOntact aims to evaluate Gene Ontology (GO) enrichments for sets of
      noncoding elements. To do this, GOntact uses chromatin contacts determined
      with the promoter capture Hi-C (PCHi-C) approach. GO categories
      are "transmitted" to noncoding elements from the genes that they
      are contacted by in the PCHi-C data. The statistical
      significance of GO enrichments is evaluated
      with binomial tests, which compare the frequency of elements
      associated with a given GO category in a foreground set with the
      frequency of elements associated with the same GO category in a
      background set.
      By default, the coordinates of enhancer elements predicted by
      ENCODE are used as a background set. For the GOntact webserver, PCHi-C chromatin contacts
      are those compiled by Laverre <i>et al.</i>, Genome Research,
      2022. The data used in this version of the GOntact webserver are
      described <a href="GOntact_data_current.html">here</a>.
      If you want to run GOntact on a different organism or a
      different set of PCHi-C samples, please see our command-line
      version on GitLab. If you want to run GREAT using the same
      genome annotations and Gene Ontology data, you can do so
      <a href="|}great_route{|">here</a>.
    </p>|}

  let%html gontact_form request = {|
  <form id= "gontact-form" method="Post" action="/" enctype="multipart/form-data">
    |}[Tyxml.Html.Unsafe.data @@ Dream.csrf_tag request]{|
    <!-- GO domain choice -->
      Select the <b>Gene Ontology domain</b> for which you want to
      compute an enrichment:
      <br>
      <input type="radio" id="radio-biological-process" name="domain-choice" value="biological_process" required/>
      <label for="radio-biological-process">biological process</label>
      <br>
      <input type="radio" id="radio-molecular-function" name="domain-choice" value="molecular_function" />
      <label for="radio-molecular-function">molecular function</label>
      <br>
       <input type="radio" id="radio-cellular-component" name="domain-choice" value="cellular_component" />
      <label for="radio-cellular-component">cellular component</label>
      <br>

      <hr>

      <!-- genome choice -->
      Select the <b>genome</b> you want to analyze:
      <br>
      <input type="radio" id="genome_human" name="genome-choice" value="human" required />
      <label for="genome_human">human (hg38)</label>
      <br>

      <input type="radio" id="genome_mouse" name="genome-choice" value="mouse" />
      <label for="genome_mouse">mouse (mm10)</label>
      <hr>

      <!-- foreground regions -->

      <label for="foreground-file">Input the coordinates of the <b>foreground</b> regions (bed format):</label>
      <br>
      <input type="file" id="foreground-file" name="foreground-file" required>

      <!-- background regions -->

      <br>

      <label for="background-file">(Optional) Input the coordinates of
      a custom set of <b>background</b> regions (bed format):</label>
      <br>
      <input type="file" id="background-file" name="background-file">

      <!-- Chromatin contact parameters-->

      <br>
      <hr>
      Selection of <b>chromatin interactions</b>:
      <br>
    	<table style="padding:0px">
	  <!-- nb samples-->
	  <tr>
	    <td>
           <label for="min-samples"> Minimum number
	of PCHi-C samples: </label>
	    </td>
	    <td>
	      <input type="number" min="1"
		     name="min-samples"  value="1" style="width: 100px">
	    </td>
	  </tr>

	   <!-- nb samples-->
	  <tr>
	    <td>
              <label for="min-score"> Minimum CHiCAGO score: </label>
	    </td>
	    <td>
	      <input type="number" min="0"
		     name="min-score"  value="5" style="width: 100px">
	    </td>
	  </tr>

	  <!-- min distance-->
	  <tr>
	    <td>
	      <label for="min-dist-contacts"> Minimum distance between
	      contacted regions (bp):</label>
	    </td>
	     <td>
	      <input type="number" min="0"
		     name="min-dist-contacts" value="25000" style="width: 100px">
	     </td>
	  </tr>

	  <!-- max distance-->
	  <tr>
	    <td>
	      <label for="max-dist-contacts"> Maximum distance between
	      contacted regions (bp):</label>
	    </td>
	     <td>
	      <input type="number" min="0"
		     name="max-dist-contacts"  value="1000000" style="width: 100px">
	     </td>
	     </tr>

	  <!-- max overlap-->
	  <tr>
	    <td>
	      <label for="max-dist-element-fragment"> Maximum distance between
	      elements and contacted fragments (bp):</label>
	    </td>
	     <td>
	      <input type="number" min="0"
		     name="max-dist-element-fragment"  value="5000" style="width: 100px">
	     </td>
	     </tr>

	  <!-- basal domain -->
	  <tr>
	    <td>
	      <label for="basal-domain"> Size of the basal regulatory domain to be
              included (bp):</label>
	    </td>
	     <td>
	      <input type="number" min="0"
		     name="basal-domain"  value="0" style="width: 100px">
	     </td>
	     </tr>
	  </table>

      <!-- submit button-->
      <div><input type="submit"></div>

    </form>|}

  let render request =
    let open Tyxml.Html in
    html_page ~mode:Form ~title:"GOntact" [
      logo_header ;
      h1 [small [txt "Gene Ontology enrichments based on chromatin contacts for noncoding elements"]] ;
      hr () ;
      intro_par ;
      hr () ;
      gontact_form request ;
    ]

  type analysis_request = {
    min_score : float ;
    basal_domain : int ;
    min_dist : int ;
    max_dist : int ;
    genome : [`human | `mouse] ;
    domain : Ontology.domain ;
    margin : int ;
    min_samples : int ;
    maybe_background_bed : string option ;
    foreground_bed : string ;
  }

  let analysis_request_decode = function
    | [ "background-file", background_file ;
        "basal-domain", basal_domain ;
        "domain-choice", domain_choice ;
        "foreground-file", foreground_file ;
        "genome-choice", genome_choice ;
        "max-dist-contacts", max_dist_contacts ;
        "max-dist-element-fragment", max_dist_element_fragment ;
        "min-dist-contacts", min_dist_contacts ;
        "min-samples", min_samples ;
        "min-score", min_score ] ->
      let open Gontact.Let_syntax.Result in
      let+ min_score = Decode_form.float "min_score" min_score
      and+ basal_domain = Decode_form.int "basal_domain" basal_domain
      and+ min_dist = Decode_form.int "min_dist" min_dist_contacts
      and+ max_dist = Decode_form.int "max_dist" max_dist_contacts
      and+ genome = Decode_form.genome "genome" genome_choice
      and+ domain = Decode_form.domain "domain" domain_choice
      and+ margin = Decode_form.int "margin" max_dist_element_fragment
      and+ min_samples = Decode_form.int "min_samples" min_samples
      and+ maybe_background_bed =
        Decode_form.maybe_single_file_contents "background_file" background_file
        |> Result.map ~f:(Option.map ~f:snd)
      and+ _, foreground_bed = Decode_form.single_file_contents "foreground_file" foreground_file
      in
      { min_score ; basal_domain ; domain ; foreground_bed ; min_dist ;
        max_dist ; genome ; margin ; min_samples ; maybe_background_bed }
    | _ -> Error `Incorrect_field_list

  let analysis
      { maybe_background_bed ; domain ; foreground_bed ;
        genome ; max_dist ; min_dist ; min_samples ; min_score ; basal_domain ;
        margin } =
    let great_param = { Great.upstream = basal_domain ; downstream = basal_domain ; extend = 0 } in
    let cea_param = { Contact_enrichment_analysis.min_score ; min_dist ; max_dist ; min_samples = Some min_samples } in
    let functional_annotation, ontology = load_functional_annotation genome ~domain in
    let gene_symbols = Functional_annotation.gene_symbols functional_annotation in
    let genome_annotation = load_genome_annotation genome ~gene_symbols in
    let annotated_baits = load_annotated_baits genome ~genome_annotation ~functional_annotation in
    let contact_graph = load_contact_graph genome ~cea_param ~annotated_baits in
    let chromosome_sizes = load_chromosome_sizes genome in
    let background_bed = Option.value_or_thunk maybe_background_bed ~default:(fun () ->
        load_default_background genome
      )
    in
    let elements =
      { FGBG.foreground = foreground_bed ; background = background_bed }
      |> FGBG.map ~f:load_bed
    in
    let enriched_terms =
      if basal_domain > 0 then
        let hea =
          Hybrid_enrichment_analysis.perform
            great_param ~contact_graph ~genome_annotation ~margin
            ~chromosome_sizes ~annotated_baits ~functional_annotation
            elements
        in
        hea.enriched_terms
      else
        let cea =
          Contact_enrichment_analysis.perform
            ~margin annotated_baits functional_annotation
            contact_graph elements
        in
        cea.enriched_terms
    in
    enriched_terms, ontology
end

module GREAT_form_page = struct
  let%html intro_par = {|
    <p style="text-align:justify;">
      This page allows you to compare the results obtained with GOntact with
      the results you would obtain
      with <a href="http://great.stanford.edu/public/html/">GREAT</a>
      using the same genomic annotations and the same set of Gene
      Ontology annotations. As for GOntact, by default, the coordinates of enhancer elements predicted by
      ENCODE are used as a background set of noncoding elements.
    </p>|}

  let%html great_form request = {|
  <form id= "great-form" method="Post" action="/great" enctype="multipart/form-data">
    |}[Tyxml.Html.Unsafe.data @@ Dream.csrf_tag request]{|
    <!-- genome choice -->
    Select the <b>Gene Ontology domain</b> for which you want to compute an enrichment:
    <br>
    <input type="radio" id="radio-biological-process" name="domain-choice" value="biological_process" required/>
    <label for="radio-biological-process">biological process</label>
    <br>
    <input type="radio" id="radio-molecular-function" name="domain-choice" value="molecular_function" />
    <label for="radio-molecular-function">molecular function</label>
    <br>
    <input type="radio" id="radio-cellular-component" name="domain-choice" value="cellular_component" />
    <label for="radio-cellular-component">cellular component</label>
    <br>

    <hr>

    <!-- genome choice -->
    Select the <b>genome</b> you want to analyze:
    <br>
    <input type="radio" id="radio-genome-human" name="genome-choice" value="human" required/>
    <label for="radio-genome-human">human (hg38)</label>
    <br>
    <input type="radio" id="radio-genome-mouse" name="genome-choice" value="mouse" />
    <label for="radio-genome-mouse">mouse (mm10)</label>
    <hr>

    <!-- foreground regions -->
    <label for="foreground-file">Input the coordinates of the <b>foreground</b> regions (bed format):</label>
    <br>
    <input type="file" id="foreground-file" name="foreground-file" required>

    <!-- background regions -->
    <br>
    <label for="background-file">(Optional) Input the coordinates of
      a custom set of <b>background</b> regions (bed format):</label>
    <br>
    <input type="file" id="background-file" name="background-file">

    <!-- GREAT domain parameters-->

    <br>
    <hr>
    Definition of <b>regulatory domains</b>:
    <br>
      <table style="padding:0px">

        <!-- nb samples-->
	<tr>
	  <td>
           <label for="upstream"> Basal regulatory domain, upstream
           region (bp): </label>
	    </td>
	    <td>
	      <input type="number" min="0"
		     name="basal-upstream"  value="5000" style="width: 100px">
	    </td>
	  </tr>

	   <!-- nb samples-->
	  <tr>
	    <td>
              <label for="downstream">  Basal regulatory domain,
              downstream region (bp):</label>
	    </td>
	    <td>
	      <input type="number" min="0"
		     name="basal-downstream"  value="1000" style="width: 100px">
	    </td>
	  </tr>

	  <!-- extension-->
	  <tr>
	    <td>
	      <label for="extension"> Regulatory domain extension size (bp):</label>
	    </td>
	     <td>
	      <input type="number" min="0"
		     name="extension" value="1000000" style="width: 100px">
	     </td>
	  </tr>



	  </table>

      <!-- submit button-->
      <div style="text-align:center"><input type="submit"></div>
    </form>|}

  let render request =
    let open Tyxml.Html in
    html_page ~mode:Form ~title:"GOntact" [
      logo_header ;
      hr () ;
      intro_par ;
      hr () ;
      great_form request ;
    ]

  type analysis_request = {
    basal_downstream : int ;
    basal_upstream : int ;
    extension : int ;
    genome : [`human | `mouse] ;
    domain : Ontology.domain ;
    maybe_background_bed : string option ;
    foreground_bed : string ;
  }

  let analysis_request_decode = function
    | [ "background-file", background_file ;
        "basal-downstream", basal_downstream ;
        "basal-upstream", basal_upstream ;
        "domain-choice", domain_choice ;
        "extension", extension ;
        "foreground-file", foreground_file ;
        "genome-choice", genome_choice ] ->
      let open Gontact.Let_syntax.Result in
      let+ basal_downstream = Decode_form.int "basal_downstream" basal_downstream
      and+ basal_upstream = Decode_form.int "basal_upstream" basal_upstream
      and+ extension = Decode_form.int "extension" extension
      and+ genome = Decode_form.genome "genome" genome_choice
      and+ domain = Decode_form.domain "domain" domain_choice
      and+ maybe_background_bed =
        Decode_form.maybe_single_file_contents "background_file" background_file
        |> Result.map ~f:(Option.map ~f:snd)
      and+ _, foreground_bed = Decode_form.single_file_contents "foreground_file" foreground_file
      in
      { basal_downstream ; basal_upstream ; domain ; foreground_bed ;
        genome ; extension ; maybe_background_bed }
    | _ -> Error `Incorrect_field_list

  let analysis
      { maybe_background_bed ; domain ; foreground_bed ;
        genome ; basal_downstream ; basal_upstream ; extension } =
    let param = { Great.upstream = basal_upstream ; downstream = basal_downstream ; extend = extension } in
    let functional_annotation, ontology = load_functional_annotation genome ~domain in
    let gene_symbols = Functional_annotation.gene_symbols functional_annotation in
    let genome_annotation = load_genome_annotation genome ~gene_symbols in
    let chromosome_sizes = load_chromosome_sizes genome in
    let background_bed = Option.value_or_thunk maybe_background_bed ~default:(fun () ->
        load_default_background genome
      )
    in
    let elements =
      { FGBG.foreground = foreground_bed ; background = background_bed }
      |> FGBG.map ~f:load_bed
    in
    let gea =
      Great.enrichment_analysis param ~genome_annotation
        ~chromosome_sizes ~functional_annotation
        elements in
    let enriched_terms = gea.enriched_terms in
    enriched_terms, ontology
end

let create_analysis_request analysis req =
  let id =
    Core_unix.gettimeofday ()
    |> Core_unix.localtime
    |> Core_unix.sexp_of_tm
    |> Sexp.to_string_mach
    |> Md5.digest_string
    |> Md5.to_hex
  in
  let (t, waiter) = Lwt.wait () in
  Hashtbl.set request_table ~key:id ~data:t ;
  Lwt.async (fun () ->
      Lwt_preemptive.detach (fun () ->
          let res = analysis req in
          Lwt.wakeup waiter res ;
        ) ()
    ) ;
  id

let html_to_string html =
  Format.asprintf "%a" (Tyxml.Html.pp ()) html

let html_get_run run_id =
  let open Tyxml.Html in
  let contents = [
    logo_header ;
    div ~a:[a_id "result-table"] []
  ]
  in
  html_page ~mode:(Results { id = run_id }) ~title:"GOntact results" contents
  |> html_to_string
  |> Dream.html

let api_get_run ~mime_type ~serializer run_id =
  let response = Dream.response ~headers:["Content-Type", mime_type] in
  match Hashtbl.find request_table run_id with
  | None -> Lwt.return (response ~status:`Not_Found "")
  | Some t ->
      let%lwt enriched_terms, ontology = t in
      let gonames = Ontology.term_names ontology in
      let compare_by_fdr x y = Go_enrichment.(Float.compare x.fdr y.fdr) in
      let ordered_results = List.sort enriched_terms ~compare:compare_by_fdr in
      List.map ordered_results ~f:(fun er ->
          let go_term = Map.find_exn gonames er.id
          and enrichment = er.observed /. er.expected in
          { go_id = er.id ; go_term ; enrichment ; pval = er.pval ; fdr = er.fdr }
        )
      |> serializer
      |> response
      |> Lwt.return

let tsv_get_run = api_get_run ~mime_type:"text/tab-separated-values" ~serializer:(fun enriched_terms ->
    List.map enriched_terms ~f:(fun er ->
        sprintf "%s\t%s\t%f\t%f\t%f" er.go_id er.go_term er.enrichment er.pval er.fdr
      )
    |> List.cons "GO id\tGO term\tenrichment\tpval\tfdr"
    |> String.concat ~sep:"\n"
  )

let json_get_run = api_get_run ~mime_type:Dream.application_json ~serializer:(fun enriched_terms ->
    enriched_terms
    |> [%yojson_of: Gontact_shared.enriched_term list]
    |> Yojson.Safe.to_string
  )


module Logger = struct
  type entry = {
    time : float ;
    user_agent : string option ;
    user_info : user_info option ;
  }
  and user_info = {
    country : string option ;
    hash : string ;
  }
  [@@deriving sexp]

  let open_db () =
    let module Project_dirs = Directories.Project_dirs (struct
        let qualifier = "fr"
        let organization = "univ-lyon1"
        let application = "gontact"
      end) in

    let data_path =
      match Project_dirs.data_dir with
      | None -> failwith "can't compute data_dir path"
      | Some path -> path
    in

    let db_path = Fpath.(data_path // (v "db") |> to_string) in

    let () =
      match Sys_unix.file_exists db_path with
      | `No -> Core_unix.mkdir_p db_path ~perm:0o750;
      | `Yes -> ()
      | `Unknown -> failwithf "please check candidate db location %s" db_path ()
    in
    LevelDB.open_db db_path

  let get_country point =
    if String.is_empty point then None
    else
      let gi = Geoip.init_exn Geoip.GEOIP_MEMORY_CACHE in
      let country_name = Geoip.country_name_by_name gi point in
      Geoip.close gi ;
      country_name

  let get_user_info ip =
    let country = get_country ip in
    let hash = Stdlib.(Digest.to_hex (Digest.string ip)) in
    { country ; hash }

  let get_ip req =
    Dream.header req "X-Forwarded-For"

  let get_user_agent req =
    Dream.header req "User-Agent"

  let add_entry db req =
    let time = Core_unix.gettimeofday () in
    let key =
      Int64.bits_of_float time
      |> Int64.to_string
    in
    let maybe_ip = get_ip req in
    let user_info = Option.map ~f:get_user_info maybe_ip in
    let user_agent = get_user_agent req in
    let entry = { time ; user_agent ; user_info } in
    let data = Sexp.to_string_mach (sexp_of_entry entry) in
    LevelDB.put db key data
end

let analysis_service ~route ~param_decode ~create_analysis_request db =
  Dream.post route (fun request ->
      match%lwt Dream.multipart request with
      | `Ok form -> (
          match param_decode form with
          | Ok req -> (
              Logger.add_entry db request ;
              let id = create_analysis_request req in
              let url = sprintf "/run/%s" id in
              Dream.redirect ~status:`See_Other request url
            )
          | Error `Incorrect_field_list ->
            print_endline ([%show: (string * (string option * string) list) list] form) ;
            Dream.empty `Bad_Request
          | Error (`Msg msg) -> Dream.html ~code:500 msg
          | Error (`Param_parsing p) -> Dream.html ~code:400 (sprintf "could not parse %s" p)
        )
      | _ -> Dream.empty `Bad_Request
    )

let () =
  let db = Logger.open_db () in
  protect ~f:(fun () ->
      Dream.run
      @@ Dream.logger
      @@ Dream.memory_sessions
      @@ Dream.router [

        Dream.get "/css/**" @@ Dream.static "server/static/css" ;
        Dream.get "/img/**" @@ Dream.static "server/static/img" ;
        Dream.get "/js/**" @@ Dream.static "server/static/js" ;

        Dream.get  "/" (fun request ->
            Dream.html (html_to_string @@ Form_page.render request)
          );

        analysis_service db ~route:"/"
          ~param_decode:Form_page.analysis_request_decode
          ~create_analysis_request:(create_analysis_request Form_page.analysis);

        Dream.get great_route (fun request ->
            Dream.html (html_to_string @@ GREAT_form_page.render request)
          );

        analysis_service db ~route:great_route
          ~param_decode:GREAT_form_page.analysis_request_decode
          ~create_analysis_request:(create_analysis_request GREAT_form_page.analysis) ;

        Dream.get "/run/:run_id"  (fun request ->
            let run_id = Dream.param request "run_id" in
            let format = Dream.query request "format" in
            match format with
            | None -> html_get_run run_id
            | Some "json" -> json_get_run run_id
            | Some "tsv"  -> tsv_get_run run_id
            | Some f ->
              let msg = sprintf "Format %s is not supported" f in
              Dream.html ~status:`Bad_Request msg
          ) ;
      ]
    )
    ~finally:(fun () -> LevelDB.close db)
