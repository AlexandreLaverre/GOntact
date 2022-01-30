open Core

let basal_domain_of_tss (gi:Genomic_interval.t) ~upstream:u ~downstream:d ~chromosome_size:cs =
  let (id, chr, start_pos, strand) = (Genomic_interval.id gi, Genomic_interval.chr gi, Genomic_interval.start_pos gi,  Genomic_interval.strand gi) in  (* start_pos is tss *)
  match strand with
  | "+" | "1" ->
    let new_start = max 1 (start_pos - u) in
    let new_end = min (start_pos + d - 1) cs in
    Genomic_interval.make ~id chr new_start new_end strand
  | "-" | "-1" ->
    let new_start = max 1 (start_pos - d + 1) in
    let new_end = min (start_pos + u) cs in  
    Genomic_interval.make ~id chr new_start new_end  strand
  | _  -> invalid_arg "gene strand should be one of +, -, 1, or -1"

let rec extend_domains_right ~ordered_tss ~init ~upstream ~downstream ~extend ~chromosome_size =
  (*basal domains are sorted by start position*)
  match ordered_tss with
  | [] -> init
  | h :: t  -> 
    match t with
    | [] ->  (* last interval on the chromosome *)
      let id = Genomic_interval.id h in
      let tss = Genomic_interval.start_pos h in (*tss pos*)
      let basal_domain = basal_domain_of_tss h ~upstream ~downstream ~chromosome_size in 
      let this_end = Genomic_interval.end_pos basal_domain in 
      let extended_end = min (tss + extend) chromosome_size in
      let new_end = max this_end extended_end in 
      (id, new_end) :: init
    | n :: _ ->
      let id = Genomic_interval.id h in
      let tss = Genomic_interval.start_pos h in (*tss pos*)
      let basal_domain = basal_domain_of_tss h ~upstream ~downstream ~chromosome_size in 
      let this_end =  Genomic_interval.end_pos basal_domain in
      let next_basal_domain = basal_domain_of_tss n ~upstream ~downstream ~chromosome_size in 
      let next_start =  Genomic_interval.start_pos next_basal_domain in
      let extended_end = min (tss + extend) (next_start - 1) in
      let new_end = max this_end extended_end in 
      let acc = (id, new_end) :: init in
      extend_domains_right ~ordered_tss:t ~init:acc ~upstream ~downstream ~extend ~chromosome_size

let rec extend_domains_left  ~reverse_ordered_tss  ~init ~upstream ~downstream ~extend ~chromosome_size =
  (*basal domains are reverse sorted by start position*)
  match reverse_ordered_tss with
  | [] -> init
  | h :: t  -> 
    match t with
    | [] ->
      let id = Genomic_interval.id h in
      let tss = Genomic_interval.start_pos h in (*tss pos*)
      let basal_domain = basal_domain_of_tss h ~upstream ~downstream ~chromosome_size in 
      let this_start = Genomic_interval.start_pos basal_domain in
      let extended_start = max (tss - extend) 1 in
      let new_start = min this_start extended_start in 
      (id, new_start) :: init
    | _ ->
      let id = Genomic_interval.id h in
      let tss = Genomic_interval.start_pos h in (*tss pos*)
      let basal_domain = basal_domain_of_tss h ~upstream ~downstream ~chromosome_size in 
      let this_start =  Genomic_interval.start_pos basal_domain in
      let max_end = Option.value_exn (List.max_elt (List.map t ~f:(fun i -> let b = (basal_domain_of_tss i ~upstream ~downstream ~chromosome_size) in Genomic_interval.end_pos b)) ~compare:compare) in (*largest end value*)
      let extended_start = max (tss - extend) (max_end + 1) in
      let new_start = min this_start extended_start in 
      let acc = (id, new_start) :: init in
      extend_domains_left ~reverse_ordered_tss:t ~init:acc ~upstream ~downstream ~extend ~chromosome_size

                  
let regulatory_domains  ~chr:chr ~chr_size:cs ~genomic_annot:ga ~upstream:u ~downstream:d ~extend:e =
  let chr_set = String.Set.of_list [chr] in 
  let filtered_annot_chr = Genomic_annotation.filter_chromosomes ga chr_set in (*take only genes on only one chromosome *)
  let major_isoforms = Genomic_annotation.identify_major_isoforms filtered_annot_chr in      (*canonical isoform for each gene*)
  let major_tss = Genomic_annotation.major_tss filtered_annot_chr ~major_isoforms in         (*genomic_interval collection TSS coordinates, they are ordered*)
  let reverse_ordered_tss = Genomic_interval_collection.reverse_sort_by_coordinate major_tss in
  let domains_right = extend_domains_right ~ordered_tss:(Genomic_interval_collection.interval_list major_tss) ~init:[] ~upstream:u ~downstream:d ~extend:e ~chromosome_size:cs in
  let domains_left = extend_domains_left ~reverse_ordered_tss:(Genomic_interval_collection.interval_list reverse_ordered_tss) ~init:[] ~upstream:u ~downstream:d ~extend:e ~chromosome_size:cs in
  let domains_right_map = String.Map.of_alist_exn domains_right in
  let domains_left_map = String.Map.of_alist_exn domains_left in
  let domains_list = List.map (String.Map.keys domains_left_map) ~f:(fun id -> Genomic_interval.make ~id:id chr (String.Map.find_exn domains_left_map id) (String.Map.find_exn domains_right_map id) ".") in 
  Genomic_interval_collection.of_interval_list domains_list

