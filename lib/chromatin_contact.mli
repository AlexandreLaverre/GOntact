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

val select_min_score : t list -> min_score:float -> t list

val select_cis : t list -> t list

val select_distance : t list -> min_dist:float -> max_dist:float -> t list
