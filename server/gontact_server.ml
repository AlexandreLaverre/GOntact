open Tyxml
open Core
open Gontact

let css_link ?a href = Html.link ?a ~rel:[`Stylesheet] ~href ()
(* let ext_js_script ?(a = []) href = Html.script ~a:(Html.a_src href :: a) (Html.txt "") *)

let html_page ~title body =
  let head =
    Html.head (Html.title (Html.txt title)) [
      Html.meta ~a:[Html.a_charset "utf-8"] () ;
      Html.meta ~a:[Html.a_http_equiv "X-UA-Compatible" ; Html.a_content "IE=edge"] () ;
      Html.meta ~a:[Html.a_name "viewport" ; Html.a_content "width=device-width, initial-scale=1"] () ;
      Html.link ~rel:[`Icon] ~href:"img/favicon.ico" () ;
      css_link "css/marx.css" ;
      css_link "css/gontact.css" ;
    ]
  in
  Html.html ~a:[Html.a_lang "en"] head (Html.body [Html.main body])

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
      <a href="GREAT.html">here</a>.
    </p>|}

  let%html gontact_form request = {|
  <form method="Post" action = "/" enctype="multipart/form-data">
    |}[Tyxml.Html.Unsafe.data @@ Dream.csrf_tag request]{|
    <!-- GO domain choice -->
      Select the <b>Gene Ontology domain</b> for which you want to
      compute an enrichment:
      <br>
      <input type="radio" id="biological_process" name="domain-choice" value="biological_process" />
      <label for="biological_process">biological process</label>
      <br>
      <input type="radio" id="molecular_function" name="domain-choice" value="molecular_function" />
      <label for="molecular_function">molecular function</label>
      <br>
       <input type="radio" id="cellular_component" name="domain-choice" value="cellular_component" />
      <label for="cellular_component">cellular component</label>
      <br>

      <hr>

      <!-- genome choice -->
      Select the <b>genome</b> you want to analyze:
      <br>
      <input type="radio" id="genome_human" name="genome-choice" value="human" />
      <label for="genome_human">human (hg38)</label>
      <br>

      <input type="radio" id="genome_mouse" name="genome-choice" value="mouse" />
      <label for="genome_mouse">mouse (mm10)</label>
      <hr>

      <!-- foreground regions -->

      <label for="foreground-file">Input the coordinates of the <b>foreground</b> regions (bed format):</label>
      <br>
      <input type="file" id="foreground-file" name="foreground-file">

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
	      <input type="text"
		     name="min-samples"  value="1" style="width: 100px">
	    </td>
	  </tr>

	   <!-- nb samples-->
	  <tr>
	    <td>
              <label for="min-score"> Minimum CHiCAGO score: </label>
	    </td>
	    <td>
	      <input type="text"
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
	      <input type="text"
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
	      <input type="text"
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
	      <input type="text"
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
	      <input type="text"
		     name="basal-domain"  value="0" style="width: 100px">
	     </td>
	     </tr>

	  </table>

      <!-- submit button-->
      <div><input type="submit"></div>

    </form>|}

  let render request =
    let open Tyxml.Html in
    html_page ~title:"GOntact" [
      header [
        img ~src:"img/GOntact_logo.png" ~alt:"GOntact logo" ~a:[a_width 600] () ;
      ] ;
      h1 [small [txt "Gene Ontology enrichments based on chromatin contacts for noncoding elements"]] ;
      hr () ;
      intro_par ;
      hr () ;
      gontact_form request ;
    ]
end

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

end

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

let analysis
    ~background_file ~domain_choice ~foreground_file
    ~genome_choice ~max_dist_contacts ~max_dist_element_fragment
    ~min_dist_contacts ~min_samples ~min_score ~basal_domain =
  let open Gontact.Let_syntax.Result in
  let+ min_score = Decode_form.float "min_score" min_score
  and+ basal_domain = Decode_form.int "basal_domain" basal_domain
  and+ min_dist = Decode_form.float "min_dist" min_dist_contacts
  and+ max_dist = Decode_form.float "max_dist" max_dist_contacts
  and+ genome = Decode_form.genome "genome" genome_choice
  and+ domain = Decode_form.domain "domain" domain_choice
  and+ margin = Decode_form.int "margin" max_dist_element_fragment
  and+ min_samples = Decode_form.int "min_samples" min_samples
  and+ _, background_bed = Decode_form.single_file_contents "background_file" background_file
  and+ _, foreground_bed = Decode_form.single_file_contents "foreground_file" foreground_file
  in
  let great_param = { Great.upstream = basal_domain ; downstream = basal_domain ; extend = 0 } in
  let cea_param = { Contact_enrichment_analysis.min_score ; min_dist ; max_dist ; min_samples = Some min_samples } in
  let functional_annotation, ontology = load_functional_annotation genome ~domain in
  let gene_symbols = Functional_annotation.gene_symbols functional_annotation in
  let genome_annotation = load_genome_annotation genome ~gene_symbols in
  let annotated_baits = load_annotated_baits genome ~genome_annotation ~functional_annotation in
  let contact_graph = load_contact_graph genome ~cea_param ~annotated_baits in
  let chromosome_sizes = load_chromosome_sizes genome in
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


let table_of_enriched_terms ers ~gonames =
  let module H = Tyxml.Html in
  let compare_by_fdr x y = Go_enrichment.(Float.compare x.fdr y.fdr) in
  let ordered_results = List.sort ers ~compare:compare_by_fdr in
  let rows = List.map ordered_results ~f:(fun er ->
      let name = String.Map.find_exn gonames er.id
      and enrichment = sprintf "%.2f" (er.observed /. er.expected)
      and pval = sprintf "%g" er.pval
      and fdr = sprintf "%g" er.fdr in
      H.tr [
        H.td [ H.txt name ] ;
        H.td [ H.txt enrichment ] ;
        H.td [ H.txt pval ] ;
        H.td [ H.txt fdr ] ;
      ]
    )
  in
  let thead = H.thead [ H.tr [
      H.th [H.txt "GO term name"] ;
      H.th [H.txt "Enrichment"] ;
      H.th [H.txt "p-value"] ;
      H.th [H.txt "FDR"] ;
    ] ]
  in
  H.table ~thead rows

let generate_result_page res_or_error =
  let module H = Tyxml.Html in
  let contents = match res_or_error with
    | Error _ -> [ H.txt "error" ]
    | Ok (enriched_terms, ontology) ->
      let gonames = Ontology.term_names ontology in
      [
        table_of_enriched_terms enriched_terms ~gonames ;
      ]
  in
  html_page ~title:"GOntact results" contents


let html_to_string html =
  Format.asprintf "%a" (Tyxml.Html.pp ()) html

let () =
  Dream.run
  @@ Dream.logger
  @@ Dream.memory_sessions
  @@ Dream.router [

    Dream.get "/css/**" @@ Dream.static "server/static/css" ;
    Dream.get "/img/**" @@ Dream.static "server/static/img" ;

    Dream.get  "/" (fun request ->
        Dream.html (html_to_string @@ Form_page.render request)
      );

    Dream.post "/" (fun request ->
        match%lwt Dream.multipart request with
        | `Ok [
            "background-file", background_file ;
            "basal-domain", basal_domain ;
            "domain-choice", domain_choice ;
            "foreground-file", foreground_file ;
            "genome-choice", genome_choice ;
            "max-dist-contacts", max_dist_contacts ;
            "max-dist-element-fragment", max_dist_element_fragment ;
            "min-dist-contacts", min_dist_contacts ;
            "min-samples", min_samples ;
            "min-score", min_score
          ] ->
          let%lwt res_or_error =
            Lwt.return @@
            analysis
              ~background_file ~domain_choice ~foreground_file
              ~genome_choice ~max_dist_contacts ~max_dist_element_fragment
              ~min_dist_contacts ~min_samples ~min_score ~basal_domain
          in
          let page = generate_result_page res_or_error in
          Dream.html (html_to_string page)
        | `Ok bad_form ->
          print_endline ([%show: (string * (string option * string) list) list] bad_form) ;
          Dream.empty `Bad_Request
        | _ -> Dream.empty `Bad_Request
      );

  ]
