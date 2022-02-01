open Core
    
type t = {
  int_list : Genomic_interval.t list ;
}
[@@deriving show]

let sort_by_coordinate l =
  List.sort l ~compare:Genomic_interval.compare_int  

let reverse_sort_by_coordinate c =
  let il = c.int_list in 
  let rl =  List.rev il in
  {int_list = rl}
  
let split_bed_line line =
  match String.split line ~on:'\t' with
  | chr :: start_pos :: end_pos :: tl ->
    let (id, strand) = match tl with
      | [] -> (String.concat ~sep:":" [chr ; start_pos ; end_pos ], ".")
      | [ id ] -> (id, ".")
      | id :: t ->
        match t with
        | _ :: strand :: _ -> (id, strand)
        | _ -> (id, ".")
    in
    let start_1based = (int_of_string start_pos) + 1 in   (*bed format: 0-based, not included*)
    let end_1based = (int_of_string end_pos)  in   (*we transform coordinates to 1-based, included*)
    Genomic_interval.make ~id chr start_1based end_1based strand
  | _ -> invalid_arg "bed file must have at least three fields " 
    
let of_bed_file path =
  let l = In_channel.read_lines path in 
  let il = List.map l ~f:split_bed_line in
  let sl = sort_by_coordinate il in
  {int_list = sl}

let of_chr_size_file path =
  let split_chr_line line =
    match String.split line ~on:'\t' with
    | chr :: size :: _ ->  Genomic_interval.make ~id:chr chr 1 (int_of_string size) "."
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


let rec intersect_lists l1 l2 = (* will return a list of (id1, id2) tuples for the interecting intervals *)
  match (l1, l2) with
  | ([], _) -> []
  | (_, []) -> []
  | (h1 :: t1, h2 :: t2) -> (
      let comp = Genomic_interval.compare h1 h2 in
      match comp with
      | Smaller_chr -> intersect_lists t1 l2 (* chr1 < chr2, there cannot be any more chr1 in the 2nd list, no intersection for h1 *)
      | Larger_chr -> intersect_lists l1 t2 (* chr1 > chr2, there can be intersection for h1, but not with h2*)
      | Smaller_no_overlap -> intersect_lists t1 l2 (*same chr, end1 < start2, there cannot be intersection for h1*)
      | Larger_no_overlap -> intersect_lists l1 t2  (*same chr, start1 > end2, there cannot be intersection for anything with h2 *)
      | Smaller_overlap | Larger_overlap | Equal  -> (Genomic_interval.id h1, Genomic_interval.id h2) :: intersect_lists l1 t2 (*there is intersection, but h1 can overlap with other elements as well*)                                                       
    )
    
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
      List.iter il  ~f:(fun i -> Printf.fprintf output "%s\t%s\t%d\t%d\t%s\n" (Genomic_interval.id i) (Genomic_interval.chr i) (Genomic_interval.start_pos i) (Genomic_interval.end_pos i) (Genomic_interval.strand i)))
        
