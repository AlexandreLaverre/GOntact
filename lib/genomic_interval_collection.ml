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

let rec fold_merge current_list (current_interval:Genomic_interval.t) (next_interval:Genomic_interval.t) next_list = 
  match Genomic_interval.merge current_interval next_interval with
  | None -> (
      match next_list with
      | [] -> (List.append current_list [current_interval; next_interval])
      | h :: t -> fold_merge (List.append current_list [current_interval]) next_interval h t
    )
  | Some merged_int -> (
      match next_list with
      | [] -> List.append current_list [merged_int]
      | h :: t -> fold_merge current_list merged_int h t
    )
    
let merge_coordinates c =
  let l = c.int_list in
  match l with
  | first :: second :: t -> (
      let merged_list = (fold_merge [] first second t) in
      {int_list = merged_list}
    )
  | _ -> {int_list = l}

