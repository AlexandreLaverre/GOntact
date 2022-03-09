open Core

type enrichment_result 

val foreground_vs_background_binom_test : go_frequencies_foreground:(int String.Map.t) -> go_frequencies_background:(int String.Map.t) -> enrichment_result list 

val write_output : enrichment_result list -> string -> unit


