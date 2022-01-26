type t = {
  id : string ;
  chr : string ;
  start_pos : int ;
  end_pos : int ;
}

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


