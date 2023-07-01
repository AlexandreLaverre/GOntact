open Tyxml

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
  <form method="Post" action = "/">
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
	      foreground regions and contacted fragments (bp):</label>
	    </td>
	     <td>
	      <input type="text"
		     name="max-dist-element-fragment"  value="5000" style="width: 100px">
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

let analysis
    ~background_file ~domain_choice ~foreground_file
    ~genome_choice ~max_dist_contacts ~max_dist_element_fragment
    ~min_dist_contacts ~min_samples ~min_score =
  (* let _background_file = match background_file with *)
  (*   | "" -> None *)
  (*   | _ -> Some background_file *)
  (* in *)
  [ background_file ;
    domain_choice ;
    foreground_file ;
    genome_choice ;
    max_dist_contacts ;
    max_dist_element_fragment ;
    min_dist_contacts ;
    min_samples ; min_score ]
  |> Lwt.return

let generate_result_page xs =
  let open Tyxml.Html in
  html_page ~title:"GOntact results" [
    ul (List.map (fun x -> li [txt x]) xs)
  ]

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
        match%lwt Dream.form request with
        | `Ok [
            "background-file", background_file ;
            "domain-choice", domain_choice ;
            "foreground-file", foreground_file ;
            "genome-choice", genome_choice ;
            "max-dist-contacts", max_dist_contacts ;
            "max-dist-element-fragment", max_dist_element_fragment ;
            "min-dist-contacts", min_dist_contacts ;
            "min-samples", min_samples ;
            "min-score", min_score
          ] ->
          let%lwt res =
            analysis
              ~background_file ~domain_choice ~foreground_file
              ~genome_choice ~max_dist_contacts ~max_dist_element_fragment
              ~min_dist_contacts ~min_samples ~min_score
          in
          let page = generate_result_page res in
          Dream.html (html_to_string page)
        | _ ->
          Dream.empty `Bad_Request);

  ]
