open Base

type mode = Form | Results of { id : string }
[@@deriving sexp]

type request_status =
  | In_progress
  | Completed of enriched_term list

and enriched_term = {
  go_id : string ;
  go_term : string ;
  enrichment : float ;
  pval : float ;
  fdr : float ;
}
[@@deriving sexp, yojson]
