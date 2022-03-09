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
  
let write_output results path =
    Out_channel.with_file path ~append:false ~f:(fun output ->
      Out_channel.output_string output "GOID\tNbElementsForeground\tTotalForeground\tObserved\tNbElementsBackground\tTotalBackground\tExpected\tPValue\tFDR\n" ; 
      List.iter results  ~f:(fun res -> Printf.fprintf output "%s\t%d\t%d\t%f\t%d\t%d\t%f\t%f\t%f\n" res.id res.count_foreground res.total_foreground res.observed res.count_background res.total_background res.expected res.pval res.fdr))

