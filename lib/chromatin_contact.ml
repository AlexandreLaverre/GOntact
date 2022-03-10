open Core

type t = {
  bait_chr : string ; 
  bait_start : int ;
  bait_end : int ;
  bait_name : string ;
  otherEnd_chr : string ;
  otherEnd_start : int ;
  otherEnd_end : int ; 
  otherEnd_name : string ;
  n_reads : int ;
  score : float ;
}

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
             
let of_ibed_file path ~strip_chr =
  let l = In_channel.read_lines path in 
  List.filter_map l ~f:(split_ibed_line ~strip_chr) 
  
let select_min_score l ~min_score =
  List.filter l ~f:Float.(fun cc -> cc.score >= min_score)

let select_cis l =
  List.filter l ~f:(fun cc -> String.equal cc.bait_chr cc.otherEnd_chr)

let compute_distance cc =
  match String.equal cc.bait_chr cc.otherEnd_chr with
  | true -> (
      let midpos1 = ((float_of_int cc.bait_start) +. (float_of_int cc.bait_end)) /. 2. in
      let midpos2 = ((float_of_int cc.otherEnd_start) +. (float_of_int cc.otherEnd_end)) /. 2. in
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
  let (chr, start_pos, end_pos) = (i.otherEnd_chr, i.otherEnd_start, i.otherEnd_end) in
  Genomic_interval.make chr start_pos end_pos Unstranded

let get_id_frag i =
  let (chr, start_pos, end_pos) = (i.otherEnd_chr, i.otherEnd_start, i.otherEnd_end) in
  Printf.sprintf "%s:%d-%d" chr start_pos end_pos

let get_id_bait i =
  let (chr, start_pos, end_pos) = (i.bait_chr, i.bait_start, i.bait_end) in
  Printf.sprintf "%s:%d-%d" chr start_pos end_pos

let compare i1 i2 =
  let id1 = Printf.sprintf "%s %s" (get_id_bait i1) (get_id_frag i1) in
  let id2 = Printf.sprintf "%s %s" (get_id_bait i2) (get_id_frag i2) in
  String.compare id1 id2
    
let select_unbaited l ~bait_collection = 
  let contacted_frag_list = List.map l ~f:contacted_fragment in
  let contacted_frag_collection = Genomic_interval_collection.of_interval_list contacted_frag_list in
  let intersection = Genomic_interval_collection.intersect contacted_frag_collection bait_collection in (* String.Map -> string list *)
  List.filter l ~f:(fun x -> not (String.Map.mem intersection (get_id_frag x)))
    
let go_annotate_baits ~bait_collection ~genome_annotation ~max_dist ~functional_annot =
  let tss_intervals = Genomic_annotation.all_tss_intervals genome_annotation max_dist in 
  let intersection = Genomic_interval_collection.intersect bait_collection tss_intervals in (*String.Map - key = bait ids ; values = list of gene_id:gene_symbol mixed id*)
  let get_terms_symbol id =
    match String.split id ~on:':' with 
    | [ _ ; symbol ] -> (
      let terms =  Functional_annotation.extract_terms functional_annot (`Symbol symbol) in
      match terms with
      | Some l -> l
      | None -> []
    )
    | _ -> []
  in
  let go_annot = String.Map.map intersection ~f:(fun l -> List.concat_map l ~f:get_terms_symbol) in
  String.Map.map go_annot ~f:(fun l -> List.dedup_and_sort ~compare:String.compare l)

let contacted_fragment_collection ~contacts =
  let all_fragments = List.map contacts ~f:contacted_fragment in 
  let unique_fragments = List.dedup_and_sort ~compare:Genomic_interval.compare_intervals all_fragments in
  Genomic_interval_collection.of_interval_list unique_fragments

let fragment_to_baits ~contacts =
  let tuples = List.map contacts ~f:(fun c -> (get_id_frag c, get_id_bait c)) in
  let sm = String.Map.of_alist_multi tuples in
  String.Map.map sm ~f:(fun l -> List.dedup_and_sort ~compare:String.compare l) (* unique list of contacted baits*)

let go_frequencies ~(element_coordinates:Genomic_interval_collection.t) ~(fragments:Genomic_interval_collection.t) ~fragment_to_baits ~annotated_baits =
  let intersection = Genomic_interval_collection.intersect element_coordinates fragments in
  let nb_elements = List.length (String.Map.keys intersection) in  
  let elbaits = String.Map.map intersection ~f:(fun l -> List.dedup_and_sort ~compare:String.compare (List.join (List.filter_map l ~f:(fun frag -> String.Map.find fragment_to_baits frag)))) in  
  let elgo = String.Map.map elbaits ~f:(fun l -> List.dedup_and_sort ~compare:String.compare (List.join (List.filter_map l ~f:(fun bait -> String.Map.find annotated_baits bait)))) in
  let elgolist = String.Map.to_alist elgo in
  let tuples = List.join (List.map elgolist ~f:(fun (id, l) -> List.map l ~f:(fun g -> (g, id)))) in
  let gomap = String.Map.of_alist_multi tuples in
  let counts = String.Map.map gomap ~f:(fun l -> List.length l) in
  let counts_with_total = String.Map.add_exn counts ~key:"total" ~data:nb_elements in
  counts_with_total


