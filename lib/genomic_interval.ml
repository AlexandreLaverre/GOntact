open Core
    
type t = {
  id : string ;
  chr : string ;
  start_pos : int ;
  end_pos : int ;
  strand : string ; 
}
[@@deriving show]

let make ?id chr start_pos end_pos strand =
  if start_pos > end_pos then  invalid_arg "end pos should be larger than start pos " ;
  match strand with
  | "1" | "-1" | "+" | "-" | "." -> (
      let id = 
        match id with
        | Some "" | None -> Printf.sprintf "%s:%d:%d" chr start_pos end_pos
        | Some i -> i
      in
      {id ; chr ; start_pos ; end_pos; strand}
    )
  | _ -> invalid_arg "strand should be one of 1, -1, +, -, ."
    
let intersect i j s  =
  match (String.equal i.chr j.chr) with
  | false -> false
  | true -> 
    let max_pos = max i.start_pos j.start_pos in
    let min_pos = min i.end_pos j.end_pos in
    match s with
    | `Stranded -> (String.equal i.strand j.strand) && (max_pos <= min_pos)
    | `Unstranded -> max_pos <= min_pos

let merge i j = 
  match (String.equal i.chr j.chr) with
  | false -> None
  | true ->
    match (String.equal i.strand j.strand) with
    | false -> None
    | true -> 
      let max_pos = max i.start_pos j.start_pos in
      let min_pos = min i.end_pos j.end_pos in
      if max_pos <= (min_pos + 1) then
        (
          let chr = i.chr in
          let start_pos = min i.start_pos j.start_pos in
          let end_pos = max i.end_pos j.end_pos in
          let strand = i.strand in
          let id = String.concat ~sep:":" [i.chr ; string_of_int start_pos ; string_of_int end_pos ] in 
          Some {id ; chr ; start_pos ; end_pos; strand}
        )
      else None
    
(* we will sort intervals by start position only - intervals may be overlapping*)
(* strand does not intervene in comparison! *)
let compare i j =
  let cmp1 = String.compare i.chr j.chr in
  if cmp1 < 0 then -3 else ( 
    if cmp1 > 0 then 3 else (
      if i.end_pos < j.start_pos then -2 else ( (*smaller, there cannot be intersection*)
        if i.start_pos < j.start_pos then -1 else ( (* smaller, but there is intersection*)
          if i.start_pos = j.start_pos then 0 else (
            if i.start_pos <= j.end_pos then 1 else ( (*larger, but there is intersection *)
              2 )))))) (*larger, no intersection possible *)

let chr t = t.chr

let id t = t.id

let start_pos t = t.start_pos

let end_pos t = t.end_pos

let strand t = t.strand


