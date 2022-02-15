type t = {
  chr1 : string ; 
  start1 : int ;
  end1 : int ;
  chr2 : string ;
  start2 : int ; 
  end2 : int ;
  n_reads : int ;
  score : float ;
}

val of_ibed_file : string -> strip_chr:bool -> t list 
