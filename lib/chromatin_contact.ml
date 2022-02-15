open Core

type t = {
  chr1 : string ; 
  start1 : int ;
  end1 : int ;
  chr2 : string ;
  start2 : int ; 
  end2 : int ;
  n_reads : int ;
  score : float ;
}

let split_ibed_line line ~strip_chr =
  if String.is_prefix line ~prefix:"bait_chr" then None (*we ignore header*)
  else
    match String.split line ~on:'\t' with
    | o_chr1 :: str_start1 :: str_end1 :: o_chr2 :: str_start2 :: str_end2 :: str_n_reads :: str_score :: _  ->
      let chr1 = (if strip_chr && (String.is_prefix o_chr1 ~prefix:"chr") then String.drop_prefix o_chr1 3 else o_chr1) in
      let chr2 = (if strip_chr && (String.is_prefix o_chr2 ~prefix:"chr") then String.drop_prefix o_chr2 3 else o_chr2) in
      let start1 = int_of_string str_start1 in
      let end1 = int_of_string str_end1 in
      let start2 = int_of_string str_start2 in
      let end2 = int_of_string str_end2 in
      let n_reads = int_of_string str_n_reads in
      let score = float_of_string str_score in
      Some {chr1 ; start1; end1 ; chr2 ; start2; end2; n_reads ; score }
    | _ -> invalid_arg "ibed file must have at least 8 tab-separated columns: chr1, start1, end1, chr2, start2, end2, n_reads,  score" 
             
let of_ibed_file path ~strip_chr =
  let l = In_channel.read_lines path in 
  List.filter_map l ~f:(split_ibed_line ~strip_chr) 
  
