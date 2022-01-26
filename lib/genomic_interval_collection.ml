open Core
    
type t = {
  int_list : Genomic_interval.t list ;
}

(*transform everything to 1-based*)

let split_bed_line line =
  match String.split line ~on:'\t' with
  | chr :: start_pos :: end_pos :: tl ->
    let id = match tl with
      | [] -> String.concat ~sep:":" [chr ; start_pos ; end_pos ]
      | id :: _ -> id
    in
    Genomic_interval.make ~id chr (int_of_string start_pos) (int_of_string end_pos) 
  | _ -> invalid_arg "bed file must have at least three fields " 
           
    
let from_bed_file path =
  let l = In_channel.read_lines path in 
  let il = List.map l ~f:split_bed_line in 
  {int_list = il}
  
