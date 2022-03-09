open Core
    
type enrichment_result = {
  observed : float ;
  expected : float ;
  count_foreground : int ;
  total_foreground : int ;
  count_background : int ;
  total_background : int ;
  pval : float ;
}

let enrichment_binom_test ~count_foreground ~total_foreground ~count_background ~total_background =
  let p = (float_of_int count_background) /. (float_of_int total_background) in
  let pval = Stats.binom_test ~n_successes:count_foreground ~n_trials:total_foreground ~p in  
  let observed = (float_of_int count_foreground) /. (float_of_int total_foreground) in
  {observed ; expected = p ; count_foreground ; total_foreground ; count_background ; total_background ; pval}


let foreground_vs_background_binom_test ~go_frequencies_foreground ~go_frequencies_background =
  let total_bg = String.Map.find_exn go_frequencies_background "total" in
  let total_fg = String.Map.find_exn go_frequencies_foreground "total" in
  let ids_bg = String.Map.keys go_frequencies_background in
  let go_ids_bg = String.Set.of_list (List.filter ids_bg ~f:(fun id -> not (String.equal id "total"))) in
  let ids_fg = String.Map.keys go_frequencies_background in
  let go_ids_fg = String.Set.of_list (List.filter ids_fg ~f:(fun id -> not (String.equal id "total"))) in
  let shared_go_ids = String.Set.to_list (String.Set.inter go_ids_bg go_ids_fg) in
  let test_one_category id = (*by construction id has to be present for both fg and bg set *)
    let count_fg= String.Map.find_exn go_frequencies_foreground id in
    let count_bg = String.Map.find_exn go_frequencies_background id in
    enrichment_binom_test ~count_foreground:count_fg ~total_foreground:total_fg ~count_background:count_bg ~total_background:total_bg
  in
  let all_tests = List.map shared_go_ids ~f:(fun id -> (id, test_one_category id)) in
  String.Map.of_alist_exn all_tests

let write_output tests path =
    Out_channel.with_file path ~append:false ~f:(fun output ->
      Out_channel.output_string output "GOID\tNbElementsForeground\tTotalForeground\tObserved\tNbElementsBackground\tTotalBackground\tExpected\tPValue\n" ; 
      String.Map.iteri tests  ~f:(fun ~key ~data -> Printf.fprintf output "%s\t%d\t%d\t%f\t%d\t%d\t%f\t%f\n" key data.count_foreground data.total_foreground data.observed data.count_background data.total_background data.expected data.pval))

