open Core
    
type t = {
  id : string ;
  chr : string ;
  start_pos : int ;
  end_pos : int ;
}
[@@deriving show]

let make ?id chr start_pos end_pos =
  if start_pos > end_pos then  invalid_arg "end pos should be larger than start pos " ;
  let id = 
    match id with
    | Some "" | None -> Printf.sprintf "%s:%d:%d" chr start_pos end_pos
    | Some i -> i
  in
  {id ; chr ; start_pos ; end_pos}

let intersect i j =
  match (String.equal i.chr j.chr) with
  | false -> false
  | true -> 
    let max_pos = max i.start_pos j.start_pos in
    let min_pos = min i.end_pos j.end_pos in
    max_pos <= min_pos

let merge i j = 
  match (String.equal i.chr j.chr) with
  | false -> None
  | true -> 
    let max_pos = max i.start_pos j.start_pos in
    let min_pos = min i.end_pos j.end_pos in
    if max_pos <= (min_pos + 1) then
      (
        let chr = i.chr in
        let start_pos = min i.start_pos j.start_pos in
        let end_pos = max i.end_pos j.end_pos in
        let id = String.concat ~sep:":" [i.chr ; string_of_int start_pos ; string_of_int end_pos ] in 
        Some {id ; chr ; start_pos ; end_pos}
      )
    else None
    
(* we will sort intervals by start position only - intervals may be overlapping*)
let compare i j =
  if (not (String.equal i.chr j.chr)) then String.compare i.chr j.chr
  else compare i.start_pos j.start_pos

let chr t = t.chr
