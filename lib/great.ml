open Core

type regulatory_domain = {
  gene_id : string ;
  gene_symbol : string ;
  chr : string ;
  tss_pos : int ;
  start_pos : int ;
  end_pos : int ;
}

let basal_domain_of_tss (tss:Genomic_interval.t) ~(genomic_annotation:Genomic_annotation.t)  ~upstream:u ~downstream:d ~chromosome_size:cs =
  let (id, chr, tss_pos, strand) = (Genomic_interval.id tss, Genomic_interval.chr tss, Genomic_interval.start_pos tss,  Genomic_interval.strand tss) in  (* start_pos is tss *)
  let gene_symbol = Genomic_annotation.gene_symbol_exn genomic_annotation id in 
  match strand with
  | Genomic_interval.Forward ->
    let new_start = max 1 (tss_pos - u) in
    let new_end = min (tss_pos + d - 1) cs in
    { gene_id = id; gene_symbol; chr ; tss_pos ; start_pos = new_start ; end_pos = new_end }
  | Genomic_interval.Reverse ->
    let new_start = max 1 (tss_pos - d + 1) in
    let new_end = min (tss_pos + u) cs in  
    { gene_id = id; gene_symbol; chr ; tss_pos ; start_pos = new_start ; end_pos = new_end }
  | Genomic_interval.Unstranded -> invalid_arg "this gene is unstranded!"
    
let extend_one_domain (d:Genomic_interval.t) ~(genomic_annotation:Genomic_annotation.t) ~left_boundary ~right_boundary ~extend ~upstream ~downstream ~chromosome_size  =
  (*d is a TSS domain domain*)
  let tss = Genomic_interval.start_pos d in
  let chr = Genomic_interval.chr d in
  let id = Genomic_interval.id d in
  let gene_symbol = Genomic_annotation.gene_symbol_exn genomic_annotation id in 
  let basal_domain = basal_domain_of_tss d ~genomic_annotation ~upstream ~downstream ~chromosome_size in  
  let current_start = basal_domain.start_pos in 
  let current_end = basal_domain.end_pos  in
  let new_start = min current_start (max (tss - extend) (left_boundary + 1)) in
  let new_end = max current_end (min (tss + extend - 1) (right_boundary - 1)) in
  { gene_id = id; gene_symbol; chr ; tss_pos = tss ; start_pos = new_start ; end_pos = new_end }

let extend_domains ~genomic_annotation ~ordered_tss ~extend ~upstream ~downstream ~chromosome_size =
  let extend_one_domain = extend_one_domain ~genomic_annotation ~extend ~upstream ~downstream ~chromosome_size in
  let basal_domain_of_tss = basal_domain_of_tss ~genomic_annotation ~upstream ~downstream ~chromosome_size in
  let rec rec_extend tss_list ~acc ~previous =
    match tss_list with
    | [] -> assert false
    | [ h ]  ->  (extend_one_domain h ~left_boundary:(previous.end_pos) ~right_boundary:(chromosome_size + 1)) :: acc  
    | i1 :: (i2 :: _ as t) ->
      let b1 = basal_domain_of_tss i1 in 
      let b2 = basal_domain_of_tss i2 in 
      let ext = extend_one_domain i1 ~left_boundary:(previous.end_pos) ~right_boundary:(b2.start_pos) in 
      rec_extend t ~acc:(ext :: acc) ~previous:b1 
  in
  match ordered_tss with
  | [] -> []
  | [ d ] -> [ extend_one_domain d ~left_boundary:1 ~right_boundary:(chromosome_size+1) ]
  | i1 :: (i2 :: _ as t) ->
    let b1 = basal_domain_of_tss i1 in
    let b2 = basal_domain_of_tss i2 in
    extend_one_domain i1 ~left_boundary:0 ~right_boundary:(b2.start_pos) :: (rec_extend t ~acc:[] ~previous:b1)

let basal_plus_extension_domains  ~chr:chr ~chromosome_size ~genomic_annotation:ga ~upstream ~downstream ~extend =
  let chr_set = String.Set.singleton chr in 
  let filtered_annot_chr = Genomic_annotation.filter_chromosomes ga chr_set in (*take only genes on only one chromosome *)
  let major_isoforms = Genomic_annotation.identify_major_isoforms filtered_annot_chr in      (*canonical isoform for each gene*)
  let major_tss = Genomic_annotation.major_isoform_tss filtered_annot_chr ~major_isoforms in         (*genomic_interval collection TSS coordinates, they are ordered*)
  let domains_list = extend_domains ~genomic_annotation:ga ~ordered_tss:(Genomic_interval_collection.interval_list major_tss) ~extend ~upstream ~downstream ~chromosome_size in
  domains_list
    
