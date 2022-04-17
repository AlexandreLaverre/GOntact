let chrono label f x =
  let t1 = Unix.time () in
  let y = f x in
  let t2 = Unix.time () in
  let elapsed = int_of_float (t2 -. t1) in
  Printf.printf "%d seconds elapsed while running '%s'.\n%!" elapsed label;
  y


