let chrono label f x =
  let t1 = Unix.time () in
  let y = f x in
  let t2 = Unix.time () in
  let elapsed = int_of_float (t2 -. t1) in
  Logs.info (fun m -> m "%d seconds elapsed while running '%s'%!" elapsed label) ;
  y

let time () = (Unix.times ()).Unix.tms_utime
let tic, tac =
  let t = ref (time ()) in
  (fun () -> t := time ()),
  (fun () ->
     let t' = time () in
     Logs.info (fun m -> m "%.1fs elapsed" (t' -. !t)) ;
     t := t')
