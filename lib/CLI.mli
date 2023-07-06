type parlist = {
  mode : mode ;
  functional_annot : string ;
  obo_path : string ;
  domain : string ;
  gene_info : string ;
  fg_path : string ;
  bg_path : string ;
  chr_sizes : string ;
  upstream : int ;
  downstream : int ;
  extend : int ;
  bait_coords : string ;
  ibed_files : string list ;
  max_dist_bait_TSS : int ;
  max_dist_element_fragment : int ;
  min_dist_contacts : int ;
  max_dist_contacts : int ;
  min_score : float ;
  output_dir : string ;
  output_prefix : string ;
  write_elements_foreground : bool ;
  write_elements_background : bool ;
}
and mode = GREAT | Contacts | Hybrid
[@@deriving sexp]

val great_mode :
  parlist ->
  chromosome_sizes:Genomic_interval_collection.t ->
  gonames:string Core.String.Map.t ->
  filtered_annot:Genomic_annotation.t ->
  foreground:Genomic_interval_collection.t ->
  background:Genomic_interval_collection.t ->
  propagated_fa:Functional_annotation.t ->
  unit

val contacts_mode :
  parlist ->
  gonames:string Core.String.Map.t ->
  filtered_annot:Genomic_annotation.t ->
  foreground:Genomic_interval_collection.t ->
  background:Genomic_interval_collection.t ->
  propagated_fa:Functional_annotation.t ->
  unit

val hybrid_mode :
  parlist ->
  chromosome_sizes:Genomic_interval_collection.t ->
  filtered_annot:Genomic_annotation.t ->
  foreground:Genomic_interval_collection.t ->
  background:Genomic_interval_collection.t ->
  propagated_fa:Functional_annotation.t ->
  gonames:string Core.String.Map.t ->
  unit

val command : unit Cmdliner.Cmd.t
