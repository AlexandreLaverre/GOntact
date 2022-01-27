open Core
    
type t = {
  int_list : Genomic_interval.t list ;
}
[@@deriving show]

let sort_by_coordinate l =
  List.sort l ~compare:Genomic_interval.compare  

let split_bed_line line =
  match String.split line ~on:'\t' with
  | chr :: start_pos :: end_pos :: tl ->
    let id = match tl with
      | [] -> String.concat ~sep:":" [chr ; start_pos ; end_pos ]
      | id :: _ -> id
    in
    let start_1based = (int_of_string start_pos) + 1 in   (*bed format: 0-based, not included*)
    let end_1based = (int_of_string end_pos)  in   (*we transform coordinates to 1-based, included*)
    Genomic_interval.make ~id chr start_1based end_1based
  | _ -> invalid_arg "bed file must have at least three fields " 
    
let from_bed_file path =
  let l = In_channel.read_lines path in 
  let il = List.map l ~f:split_bed_line in
  let sl = sort_by_coordinate il in
  {int_list = sl}

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

