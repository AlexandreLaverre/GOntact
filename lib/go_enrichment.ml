(*open Core*)

open OCamlR_stats

type enrichment_result = {
  observed : float ;
  expected : float ;
  binom_pval : float ;
}

(*

let binom_test ~count_foreground ~total_count_foreground ~count_background ~total_count_background =
  
  *)
(*
let foreground_vs_background ~go_frequencies_foreground  ~total_count_foreground ~go_frequencies_background ~total_count_background  =
  *)
