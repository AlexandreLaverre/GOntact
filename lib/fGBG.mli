type 'a t = {
  foreground : 'a ;
  background : 'a ;
}

val map : 'a t -> f:('a -> 'b) -> 'b t
