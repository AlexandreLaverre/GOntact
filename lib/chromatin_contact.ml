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
  
let select_min_score l ~min_score =
  List.filter l ~f:Float.(fun cc -> cc.score >= min_score)

let select_cis l =
  List.filter l ~f:(fun cc -> String.equal cc.chr1 cc.chr2)

let compute_distance cc =
  match String.equal cc.chr1 cc.chr2 with
  | true -> (
      let midpos1 = ((float_of_int cc.start1) +. (float_of_int cc.end1)) /. 2. in
      let midpos2 = ((float_of_int cc.start2) +. (float_of_int cc.end2)) /. 2. in
      let dist = Float.abs (midpos1 -. midpos2) in
      Some dist
    )
  | false -> None

let select_distance l ~min_dist ~max_dist =
  let verify_distance cc =
    let d = compute_distance cc in
    match d with
    | None -> false
    | Some dist -> Float.(dist >= min_dist && dist <= max_dist) 
  in      
  List.filter l ~f:verify_distance

let contacted_fragment i =
  let (chr, start_pos, end_pos) = (i.chr2, i.start2, i.end2) in
  Genomic_interval.make chr start_pos end_pos Unstranded

(*
let get_id_bait i =
  let (chr, start_pos, end_pos) = (i.chr1, i.start1, i.end1) in
  Printf.sprintf "%s:%d:%d" chr start_pos end_pos
*)

let get_id_frag i =
  let (chr, start_pos, end_pos) = (i.chr2, i.start2, i.end2) in
  Printf.sprintf "%s:%d:%d" chr start_pos end_pos
  
let select_unbaited l ~bait_collection = 
  let contacted_frag_list = List.map l ~f:contacted_fragment in
  let contacted_frag_collection = Genomic_interval_collection.of_interval_list contacted_frag_list in
  let intersection = Genomic_interval_collection.intersect contacted_frag_collection bait_collection in (* String.Map -> string list *)
  List.filter l ~f:(fun x -> not (String.Map.mem intersection (get_id_frag x)))
    
