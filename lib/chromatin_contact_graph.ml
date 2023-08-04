open Core
open Chromatin_contact

type t = Chromatin_contact.t list

let split_ibed_line line ~strip_chr =
  if String.is_prefix line ~prefix:"bait_chr" then None (*we ignore header*)
  else
    match String.split line ~on:'\t' with
    | [ o_chr1 ; str_start1 ; str_end1 ; name1 ;  o_chr2 ;  str_start2 ; str_end2 ; name2 ; str_n_reads ; str_score ]  ->
      let chr1 = (if strip_chr && (String.is_prefix o_chr1 ~prefix:"chr") then String.drop_prefix o_chr1 3 else o_chr1) in
      let chr2 = (if strip_chr && (String.is_prefix o_chr2 ~prefix:"chr") then String.drop_prefix o_chr2 3 else o_chr2) in
      let start1 = int_of_string str_start1 in
      let end1 = int_of_string str_end1 in
      let start2 = int_of_string str_start2 in
      let end2 = int_of_string str_end2 in
      let n_reads = int_of_string str_n_reads in
      let score = float_of_string str_score in
      Some {bait_chr = chr1 ; bait_start = start1; bait_end = end1 ; bait_name = name1 ; otherEnd_chr = chr2 ; otherEnd_start = start2; otherEnd_end = end2; otherEnd_name = name2 ; n_reads ; score }
    | _ -> invalid_arg "ibed file must have 10 tab-separated columns: bait_chr, bait_start, bait_end, bait_name, otherEnd_chr, otherEnd_start, otherEnd_end, otherEnd_name, n_reads,  score"

let split_ibed_line_filtered line ~strip_chr ~min_dist ~max_dist ~min_score ~bait_map  =
  if String.is_prefix line ~prefix:"bait_chr" then None (*we ignore header*)
  else
    match String.split line ~on:'\t' with
    | [ o_chr1 ; str_start1 ; str_end1 ; name1 ;  o_chr2 ;  str_start2 ; str_end2 ; name2 ; str_n_reads ; str_score ]  ->
      let chr1 = (if strip_chr && (String.is_prefix o_chr1 ~prefix:"chr") then String.drop_prefix o_chr1 3 else o_chr1) in
      let chr2 = (if strip_chr && (String.is_prefix o_chr2 ~prefix:"chr") then String.drop_prefix o_chr2 3 else o_chr2) in
      let start1 = int_of_string str_start1 in
      let end1 = int_of_string str_end1 in
      let start2 = int_of_string str_start2 in
      let end2 = int_of_string str_end2 in
      let n_reads = int_of_string str_n_reads in
      let score = float_of_string str_score in
      if String.equal chr1 chr2 then
        if Float.(score >= min_score) then
          let midpos1 = ((float_of_int start1) +. (float_of_int end1)) /. 2. in
          let midpos2 = ((float_of_int start2) +. (float_of_int end2)) /. 2. in
          let dist = Float.abs (midpos1 -. midpos2) in
          if Float.(dist >= min_dist && dist <= max_dist) then
            let id_frag = Printf.sprintf "%s:%d-%d" chr2 start2 end2 in
            if String.Set.mem bait_map id_frag then None
            else Some {bait_chr = chr1 ; bait_start = start1; bait_end = end1 ; bait_name = name1 ; otherEnd_chr = chr2 ; otherEnd_start = start2; otherEnd_end = end2; otherEnd_name = name2 ; n_reads ; score }
          else None
        else None
      else None
    | _ -> invalid_arg "ibed file must have 10 tab-separated columns: bait_chr, bait_start, bait_end, bait_name, otherEnd_chr, otherEnd_start, otherEnd_end, otherEnd_name, n_reads,  score"

let of_ibed_file path ~strip_chr =
  let l = In_channel.read_lines path in
  List.filter_map l ~f:(split_ibed_line ~strip_chr)

let of_ibed_file_filtered path ~strip_chr ~min_dist ~max_dist ~min_score ~bait_map =
  let l = In_channel.read_lines path in
  List.filter_map l ~f:(split_ibed_line_filtered ~strip_chr ~min_dist ~max_dist ~min_score ~bait_map)

let select_min_score l ~min_score =
  List.filter l ~f:Float.(fun cc -> cc.score >= min_score)

let select_cis l =
  List.filter l ~f:(fun cc -> String.equal cc.bait_chr cc.otherEnd_chr)

let select_distance l ~min_dist ~max_dist =
  let verify_distance cc =
    let d = Chromatin_contact.distance cc in
    match d with
    | None -> false
    | Some dist -> dist >= min_dist && dist <= max_dist
  in
  List.filter l ~f:verify_distance

let select_unbaited l ~bait_collection =
  let bait_ids = String.Set.of_list (List.map (Genomic_interval_collection.interval_list bait_collection) ~f:(fun i -> Genomic_interval.id i)) in
  let unbaited = List.filter l ~f:(fun x -> not (String.Set.mem bait_ids (get_id_frag x))) in
  unbaited

let remove_unannotated_baits contacts ~bait_annotation =
  List.filter contacts ~f:(fun c -> String.Map.mem bait_annotation (get_id_bait c))

let extend_fragments contacts ~margin =
  let all_fragments = List.map contacts ~f:Chromatin_contact.contacted_fragment in
  let unique_fragments = List.dedup_and_sort ~compare:Genomic_interval.compare_intervals all_fragments in
  let extended_fragments = List.map unique_fragments ~f:(fun i -> Genomic_interval.make ~id:(Genomic_interval.id i) (Genomic_interval.chr i) ((Genomic_interval.start_pos i) - margin) ((Genomic_interval.end_pos i) + margin) (Genomic_interval.strand i)) in
  Genomic_interval_collection.of_interval_list extended_fragments

let contacted_fragment_collection contacts =
  let all_fragments = List.map contacts ~f:contacted_fragment in
  let unique_fragments = List.dedup_and_sort ~compare:Genomic_interval.compare_intervals all_fragments in
  Genomic_interval_collection.of_interval_list unique_fragments

let fragment_to_baits contacts =
  let tuples = List.map contacts ~f:(fun c -> (get_id_frag c, get_id_bait c)) in
  let sm = String.Map.of_alist_multi tuples in
  String.Map.map sm ~f:(fun l -> List.dedup_and_sort ~compare:String.compare l) (* unique list of contacted baits for each fragment*)
