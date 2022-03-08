open Core

type enrichment_result 

val binom_test : count_foreground:int -> total_foreground:int -> count_background:int -> total_background:int -> enrichment_result

val foreground_vs_background_binom_test : go_frequencies_foreground:(int String.Map.t) -> go_frequencies_background:(int String.Map.t) -> (enrichment_result String.Map.t) 

val write_output : (enrichment_result String.Map.t) -> string -> unit


