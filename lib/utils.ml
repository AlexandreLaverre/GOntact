let time () = (Unix.times ()).Unix.tms_utime

let chrono label f x =
  let t1 = time () in
  let y = f x in
  let t2 = time () in
  Logs.debug (fun m -> m "%.1f seconds elapsed while running '%s'%!" (t2 -. t1) label) ;
  y

let tic, tac =
  let t = ref (time ()) in
  (fun () -> t := time ()),
  (fun () ->
     let t' = time () in
     Logs.info (fun m -> m "%.1fs elapsed" (t' -. !t)) ;
     t := t')
