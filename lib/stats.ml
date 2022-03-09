open Core

let binom_test ~n_successes ~n_trials ~p =
  1. -. Gsl.Cdf.binomial_P ~k:(n_successes - 1) ~n:n_trials ~p

let%expect_test _ =
  print_endline {|
> binom.test(4, 5, 0.2,alternative="greater")

       Exact binomial test

data:  4 and 5
number of successes = 4, number of trials = 5, p-value = 0.00672
alternative hypothesis: true probability of success is greater than 0.2
95 percent confidence interval:
 0.3425917 1.0000000
sample estimates:
probability of success
                   0.8
|} ;
  Printf.printf "pvalue = %f\n" (binom_test ~n_successes:4 ~n_trials:5 ~p:0.2) ;
  [%expect {|
    > binom.test(4, 5, 0.2,alternative="greater")

           Exact binomial test

    data:  4 and 5
    number of successes = 4, number of trials = 5, p-value = 0.00672
    alternative hypothesis: true probability of success is greater than 0.2
    95 percent confidence interval:
     0.3425917 1.0000000
    sample estimates:
    probability of success
                       0.8

    pvalue = 0.006720 |}]


let compare_named_pvalues (_, p1) (_, p2) =
  Float.compare p1 p2
                
let false_discovery_rates named_pvalues =
  let sorted_pvalues = List.sort named_pvalues ~compare:compare_named_pvalues in
  let reverse_sorted_pvalues = List.rev sorted_pvalues in
  let pval_array = Array.of_list reverse_sorted_pvalues in
  let n = Array.length pval_array in 
  let nf = float_of_int n in 
  let compute_fdr current_index current_fdr (current_id, current_pval) =
    if current_index > 0 then
      let possible_fdr = current_pval *. nf /. (float_of_int current_index) in        
      match current_fdr with
      | (_, last_fdr) :: _ -> 
        let cum_min = Float.min last_fdr possible_fdr in
        let this_fdr = Float.min cum_min 1.0 in
        (current_id, this_fdr) :: current_fdr
      | [] ->
        let this_fdr = Float.min possible_fdr 1.0 in
        [ (current_id, this_fdr) ]
    else current_fdr
  in
  Array.foldi pval_array ~init:[] ~f:(fun i current_fdr current_pval -> compute_fdr (n-i) current_fdr current_pval)
