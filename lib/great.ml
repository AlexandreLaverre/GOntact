open Core

let basal_domain_of_tss (gi:Genomic_interval.t) ~upstream:u ~downstream:d ~chromosome_size:cs =
  let (id, chr, start_pos, strand) = (Genomic_interval.id gi, Genomic_interval.chr gi, Genomic_interval.start_pos gi,  Genomic_interval.strand gi) in  (* start_pos is tss *)
  match strand with
  | Genomic_interval.Forward ->
    let new_start = max 1 (start_pos - u) in
    let new_end = min (start_pos + d - 1) cs in
    Genomic_interval.make ~id chr new_start new_end Genomic_interval.Unstranded
  | Genomic_interval.Reverse ->
    let new_start = max 1 (start_pos - d + 1) in
    let new_end = min (start_pos + u) cs in  
    Genomic_interval.make ~id chr new_start new_end Genomic_interval.Unstranded
  | Genomic_interval.Unstranded -> invalid_arg "this gene is unstranded!"
    
let extend_one_domain (d:Genomic_interval.t) ~left_boundary ~right_boundary ~extend ~upstream ~downstream ~chromosome_size  =
  (*d is a tss interval*)
  let tss = Genomic_interval.start_pos d in
  let chr = Genomic_interval.chr d in
  let id = Genomic_interval.id d in
  let basal_domain = basal_domain_of_tss d ~upstream ~downstream ~chromosome_size in  
  let current_start = Genomic_interval.start_pos basal_domain in 
  let current_end = Genomic_interval.end_pos basal_domain in
  let new_start = min current_start (max (tss - extend) (left_boundary + 1)) in
  let new_end = max current_end (min (tss + extend-1) (right_boundary - 1)) in
  Genomic_interval.make ~id  chr new_start new_end Genomic_interval.Unstranded

let extend_domains ~ordered_tss ~extend ~upstream ~downstream ~chromosome_size =
  let extend_one_domain = extend_one_domain ~extend ~upstream ~downstream ~chromosome_size in
  let basal_domain_of_tss = basal_domain_of_tss ~upstream ~downstream ~chromosome_size in
  let rec rec_extend tss_list ~acc ~previous =
    match tss_list with
    | [] -> assert false
    | [ h ]  ->  extend_one_domain h ~left_boundary:(Genomic_interval.end_pos previous) ~right_boundary:chromosome_size :: acc  
    | i1 :: (i2 :: _ as t) ->
      let b1 = basal_domain_of_tss i1 in 
      let b2 = basal_domain_of_tss i2 in 
      let ext = extend_one_domain i1 ~left_boundary:(Genomic_interval.end_pos previous) ~right_boundary:(Genomic_interval.start_pos b2) in 
      rec_extend t ~acc:(ext :: acc) ~previous:b1 
  in
  match ordered_tss with
  | [] -> []
  | [ d ] -> [ extend_one_domain d ~left_boundary:1 ~right_boundary:chromosome_size ]
  | i1 :: (i2 :: _ as t) ->
    let b1 = basal_domain_of_tss i1 in
    let b2 = basal_domain_of_tss i2 in
    extend_one_domain i1 ~left_boundary:1 ~right_boundary:(Genomic_interval.start_pos b2) :: (rec_extend t ~acc:[] ~previous:b1)

let basal_plus_extension_domains  ~chr:chr ~chromosome_size ~genomic_annot:ga ~upstream ~downstream ~extend =
  let chr_set = String.Set.singleton chr in 
  let filtered_annot_chr = Genomic_annotation.filter_chromosomes ga chr_set in (*take only genes on only one chromosome *)
  let major_isoforms = Genomic_annotation.identify_major_isoforms filtered_annot_chr in      (*canonical isoform for each gene*)
  let major_tss = Genomic_annotation.major_isoform_tss filtered_annot_chr ~major_isoforms in         (*genomic_interval collection TSS coordinates, they are ordered*)
  let domains_list = extend_domains ~ordered_tss:(Genomic_interval_collection.interval_list major_tss) ~extend ~upstream ~downstream ~chromosome_size in
  Genomic_interval_collection.of_interval_list domains_list
    
