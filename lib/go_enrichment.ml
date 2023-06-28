open Core

type enrichment_result = {
  id : string ;
  observed : float ;
  expected : float ;
  count_foreground : int ;
  total_foreground : int ;
  count_background : int ;
  total_background : int ;
  pval : float ;
  fdr : float ;
}

let compare er1 er2 =
  Float.compare er1.fdr er2.fdr

let combine_maps gomap1 gomap2 =
  let goids1 = String.Map.keys gomap1 in
  let goids2 = String.Map.keys gomap2 in
  let goids = List.append goids1 goids2 in
  let unique_goids = List.dedup_and_sort ~compare:String.compare goids in
  let gotuples = List.map unique_goids ~f:(fun goid ->
      let olist1 = String.Map.find gomap1 goid in
      let olist2 = String.Map.find gomap2 goid in
      let l = (
        match (olist1, olist2) with
        | (Some l1, Some l2) -> List.append l1 l2
        | (Some l1, _) -> l1
        | (_, Some l2) -> l2
        | (None, None) -> invalid_arg "id has to be in one of the lists!"
      ) in
      let ul = List.dedup_and_sort ~compare:String.compare l in
      (goid, ul)
    ) in
  let gomap = String.Map.of_alist_exn gotuples in
  gomap

let combine_GO_maps xs ys =
  List.map2_exn xs ys ~f:(fun (gi_x, cats_x) (gi_y, cats_y) ->
      assert Genomic_interval.(String.equal (id gi_x) (id gi_y)) ;
      gi_x, Great.GO_term_set.union cats_x cats_y
    )

let go_frequencies_legacy ~categories_by_element =
  let nb_total = String.Map.length categories_by_element in (* number of elements that have at least one GO category*)
  let counts =
    let table = String.Table.create () in
    String.Map.iter categories_by_element ~f:(fun categories ->
        List.iter categories ~f:(fun cat ->
            String.Table.update table cat ~f:(function
                | None -> 1
                | Some c -> c + 1
              )
          )
      ) ;
    String.Table.to_alist table
    |> String.Map.of_alist_exn
  in
  String.Map.add_exn counts ~key:"total" ~data:nb_total

let go_frequencies ~categories_by_element fa =
  let nb_total = List.length categories_by_element in (* number of elements that have at least one GO category*)
  let counts =
    let table = Functional_annotation.create_term_table fa 0 in
    List.iter categories_by_element ~f:(fun (_, categories) ->
        Great.GO_term_set.iter categories ~f:(fun cat ->
            let c = Ontology.PKey.Table.get table cat in
            Ontology.PKey.Table.set table cat (c + 1)
          )
      ) ;
    Ontology.PKey.Table.fold table ~init:String.Map.empty ~f:(fun ~term ~value acc ->
        if value > 0 then String.Map.add_exn acc ~key:term.id ~data:value
        else acc
      )
  in
  String.Map.add_exn counts ~key:"total" ~data:nb_total

let foreground_vs_background_binom_test ~go_frequencies_foreground ~go_frequencies_background =
  let total_bg = String.Map.find_exn go_frequencies_background "total" in
  let total_fg = String.Map.find_exn go_frequencies_foreground "total" in
  let ids_bg = String.Map.keys go_frequencies_background in
  let go_ids_bg = String.Set.of_list (List.filter ids_bg ~f:(fun id -> not (String.equal id "total"))) in
  let ids_fg = String.Map.keys go_frequencies_foreground in
  let go_ids_fg = String.Set.of_list (List.filter ids_fg ~f:(fun id -> not (String.equal id "total"))) in
  let shared_go_ids = String.Set.to_list (String.Set.inter go_ids_bg go_ids_fg) in
  let test_one_category id = (*by construction id has to be present for both fg and bg set *)
    let count_fg = String.Map.find_exn go_frequencies_foreground id in
    let count_bg = String.Map.find_exn go_frequencies_background id in
    let p = (float_of_int count_bg) /. (float_of_int total_bg) in
    let pval = Stats.binom_test ~n_successes:count_fg ~n_trials:total_fg ~p in
    (id, pval)
  in
  let all_pvalues = List.map shared_go_ids ~f:test_one_category in
  let all_pvalues_dict = String.Map.of_alist_exn all_pvalues in
  let all_fdr = Stats.false_discovery_rates all_pvalues in
  let all_fdr_dict = String.Map.of_alist_exn all_fdr in
  let construct_enrichment_object id =
    let count_fg = String.Map.find_exn go_frequencies_foreground id in
    let count_bg = String.Map.find_exn go_frequencies_background id in
    let pval = String.Map.find_exn all_pvalues_dict id in
    let fdr = String.Map.find_exn all_fdr_dict id in
    let observed = (float_of_int count_fg) /. (float_of_int total_fg) in
    let expected = (float_of_int count_bg) /. (float_of_int total_bg) in
    {id ; observed ; expected ; count_foreground = count_fg ; total_foreground = total_fg ; count_background = count_bg ; total_background = total_bg ; pval ; fdr }
  in
  List.map shared_go_ids ~f:construct_enrichment_object

let write_detailed_association gomap path =
  Out_channel.with_file path ~append:false ~f:(fun output ->
      String.Map.iteri gomap  ~f:(fun ~key ~data ->
          Printf.fprintf output "%s\t" key ;
          List.iter data ~f:(fun g -> Printf.fprintf output "%s," g) ;
          Printf.fprintf output "\n" ))

let write_output results go_names path =
  let ordered_results = List.sort results ~compare in
  Out_channel.with_file path ~append:false ~f:(fun output ->
      Out_channel.output_string output "GOID\tGOName\tNbElementsForeground\tTotalForeground\tObserved\tNbElementsBackground\tTotalBackground\tExpected\tPValue\tFDR\n" ;
      List.iter ordered_results  ~f:(fun res ->
          let name = String.Map.find_exn go_names res.id in
          Printf.fprintf output "%s\t%s\t%d\t%d\t%f\t%d\t%d\t%f\t%f\t%f\n" res.id name res.count_foreground res.total_foreground res.observed res.count_background res.total_background res.expected res.pval res.fdr))
