type 'a t = {
  foreground : 'a ;
  background : 'a ;
}

let map { foreground ; background } ~f =
  { foreground = f foreground ;
    background = f background }
