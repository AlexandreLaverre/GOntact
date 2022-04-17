open Core
    
type t = {
  int_list : Genomic_interval.t list ;
}
[@@deriving show]

type format = Base1 | Base0

let sort_by_coordinate l =
  List.sort l ~compare:Genomic_interval.compare_intervals  

let reverse_sort_by_coordinate c =
  let il = c.int_list in 
  let rl =  List.rev il in
  {int_list = rl}
  
let split_bed_line line ~strip_chr ~format =
  (*typical bed format: 0-based, not included*)
  (*if needed we transform coordinates to 1-based, included*)
  let shift = match format with
    | Base1 -> 0
    | Base0 -> 1
  in      
  match String.split line ~on:'\t' with
  | chr :: start_pos :: end_pos :: tl ->
    let (id, strand) = match tl with
      | [] -> ((Printf.sprintf "%s:%s-%s" chr start_pos end_pos), ".")
      | [ id ] -> (id, ".")
      | id :: t ->
        match t with
        | _ :: strand :: _ -> (id, strand)
        | _ -> (id, ".")
    in
    let start_1based = (int_of_string start_pos) + shift in  
    let end_1based = (int_of_string end_pos)  in   
    let new_strand = (match strand with
      | "+" -> Genomic_interval.Forward
      | "-" -> Genomic_interval.Reverse
      | "." -> Genomic_interval.Unstranded
      | _ -> invalid_arg "strand should be one of +, - or ."
      ) in
    let new_chr = (if strip_chr && (String.is_prefix chr ~prefix:"chr") then String.drop_prefix chr 3 else chr) in
    Genomic_interval.make ~id new_chr start_1based end_1based new_strand
  | _ -> invalid_arg "bed file must have at least three fields " 
    
let of_bed_file path ~strip_chr ~format =
  let l = In_channel.read_lines path in 
  let il = List.map l ~f:(split_bed_line ~strip_chr ~format) in
  let sl = sort_by_coordinate il in
  {int_list = sl}

let of_chr_size_file path ~strip_chr =
  let split_chr_line line =
    match String.split line ~on:'\t' with
    | chr :: size :: _ ->
      let new_chr = (if strip_chr && (String.is_prefix chr ~prefix:"chr") then String.drop_prefix chr 3 else chr) in
      Genomic_interval.make ~id:chr new_chr 1 (int_of_string size) Genomic_interval.Unstranded
    | _ -> invalid_arg "chr size file must have at least two fields: chr_name chr_size " 
  in
  let l = In_channel.read_lines path in
  let il = List.map l ~f:split_chr_line in
  let sl = sort_by_coordinate il in
  {int_list = sl}

let of_interval_list il = 
  let sl = sort_by_coordinate il in
  {int_list = sl}

let chr_set t = 
  let l = t.int_list in
  let chr_list = List.map l ~f:(fun x -> Genomic_interval.chr x) in
  String.Set.of_list chr_list

let fold_merge l ~init ~f = 
  let rec aux_merge acc current_interval next_list =
    match next_list with
    | h :: t -> (
        match Genomic_interval.merge current_interval h with
        | None -> aux_merge (f acc current_interval) h t 
        | Some merged_int -> aux_merge acc merged_int t
      )
    | [] -> f acc current_interval
  in
  match l with
  | h :: t ->  aux_merge init h t
  | [] -> [] 
  
let merge_coordinates c =
  let merged_list = fold_merge c.int_list ~init:[] ~f:(fun l i -> i :: l) in
  let reordered = List.rev merged_list in
  {int_list = reordered}

let intersect_lists l1 l2 =
  let chr_tuples1 = List.map l1 ~f:(fun i -> (Genomic_interval.chr i, i)) in
  let chr_map1 = String.Map.of_alist_multi chr_tuples1 in
  let chr_tuples2 = List.map l2 ~f:(fun i -> (Genomic_interval.chr i, i)) in
  let chr_map2 = String.Map.of_alist_multi chr_tuples2 in
  let chr1 = String.Set.of_list (String.Map.keys chr_map1) in
  let chr2 = String.Set.of_list (String.Map.keys chr_map2) in
  let common_chr = String.Set.to_list (String.Set.inter chr1 chr2) in
  let rec intersect_arrays a1 a2 i j startj init =
    (* Printf.printf "i %d j %d startj %d\n" i j startj ; *)
    let l1 = Array.length a1 in
    let l2 = Array.length a2 in
    if (i < l1 && j < l2) then 
      let el1 = Array.get a1 i in
      let el2 = Array.get a2 j in
      let comp = Genomic_interval.check_overlap el1 el2 in
      match comp with
      | Smaller_no_overlap -> intersect_arrays a1 a2 (i+1) startj startj init (*same chr, end1 < start2, there cannot be intersection for el1*)
      | Larger_no_overlap -> intersect_arrays a1 a2 i (j+1) (j+1) init (*same chr, start1 > end2, there cannot be intersection for anything with el2 *)
      | Smaller_overlap | Larger_overlap | Equal  ->
        (*  Printf.printf "found intersection for i %d and j %d\n" i j ; *)
        let new_init = (Genomic_interval.id el1, Genomic_interval.id el2) :: init in 
        intersect_arrays a1 a2 i (j+1) startj new_init 
      | _ -> invalid_arg "we only apply this on intervals on the same chromosome"
    else
      init
  in
  let intersection_list = List.map common_chr ~f:(fun chr ->
      let achr1 =  Array.of_list (String.Map.find_exn chr_map1 chr) in
      Array.sort achr1 ~compare:Genomic_interval.compare_intervals ;
      let achr2 =  Array.of_list (String.Map.find_exn chr_map2 chr) in
      Array.sort achr2 ~compare:Genomic_interval.compare_intervals ;
      intersect_arrays achr1 achr2 0 0 0 []
    ) in
   List.join intersection_list 
        
let intersect c1 c2 =
  let (l1, l2) = (c1.int_list, c2.int_list) in (*collections of genomic intervals are necessarily ordered *)
  let tuple_list = intersect_lists l1 l2 in
  String.Map.of_alist_multi tuple_list

let interval_list c = c.int_list

let map c ~f =
  let int_list = c.int_list in
  let l = List.map int_list ~f:f in
  of_interval_list l

let iter c ~f =
  let int_list = c.int_list in
  List.iter int_list ~f

let write_output c path ~append =
  let il = c.int_list in
  Out_channel.with_file path ~append ~f:(fun output -> 
      List.iter il  ~f:(fun i -> Printf.fprintf output "%s\t%s\t%d\t%d\t%s\n" (Genomic_interval.id i) (Genomic_interval.chr i) (Genomic_interval.start_pos i) (Genomic_interval.end_pos i) Genomic_interval.(string_of_strand (strand i))))

let length c = List.length c.int_list

let remove_duplicated_identifiers c = 
  let il = c.int_list in
  let tuples = List.map il ~f:(fun i -> (Genomic_interval.id i, i)) in
  let multi_map = String.Map.of_alist_multi tuples in
  let filtered_interval_list = String.Map.fold multi_map ~init:[] ~f:(fun ~key:_ ~data acc -> (
        match data with
        | [ int ] -> int :: acc
        | _ -> acc
      )
    ) in
  of_interval_list filtered_interval_list

let interval_map c =
  let il = c.int_list in
  let tuples = List.map il ~f:(fun i -> (Genomic_interval.id i, i)) in
  String.Map.of_alist_exn tuples (* identifiers have to be unique *)
