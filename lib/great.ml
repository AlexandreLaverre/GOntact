open Core

type regulatory_domain = {
  gene_id : string ;
  gene_symbol : string ;
  chr : string ;
  tss_pos : int ;
  start_pos : int ;
  end_pos : int ;
}

let genomic_interval_collection rdl =
  let interval_of_domain rd =
    Genomic_interval.make ~id:rd.gene_symbol rd.chr rd.start_pos rd.end_pos Genomic_interval.Unstranded
  in
  let il = List.map rdl ~f:interval_of_domain in
  Genomic_interval_collection.of_interval_list il


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

let basal_plus_extension_domains_one_chr  ~chr ~chromosome_size ~genomic_annotation:ga ~upstream ~downstream ~extend =
  let chr_set = String.Set.singleton chr in
  let filtered_annot_chr = Genomic_annotation.filter_chromosomes ga chr_set in (*take only genes on only one chromosome *)
  let major_isoforms = Genomic_annotation.identify_major_isoforms filtered_annot_chr in      (*canonical isoform for each gene*)
  let major_tss = Genomic_annotation.major_isoform_tss filtered_annot_chr ~major_isoforms in         (*genomic_interval collection TSS coordinates, they are ordered*)
  let domains_list = extend_domains ~genomic_annotation:ga ~ordered_tss:(Genomic_interval_collection.interval_list major_tss) ~extend ~upstream ~downstream ~chromosome_size in
  domains_list

let basal_plus_extension_domains ~(chromosome_sizes:Genomic_interval_collection.t) ~genomic_annotation ~upstream ~downstream ~extend =
  let cl = Genomic_interval_collection.interval_list chromosome_sizes in
  List.concat_map cl ~f:(fun i -> basal_plus_extension_domains_one_chr ~chr:(Genomic_interval.chr i) ~chromosome_size:(Genomic_interval.end_pos i) ~genomic_annotation ~upstream ~downstream ~extend)

let sorted_list_union xs ys ~compare =
  let rec loop acc xs ys =
    match xs, ys with
    | [], ys -> List.rev_append acc ys
    | xs, [] -> List.rev_append acc xs
    | hx :: tx, hy :: ty ->
      match compare hx hy with
      |  0 -> loop (hx :: acc) tx ty
      | -1 -> loop (hx :: acc) tx ys
      |  1 -> loop (hy :: acc) xs ty
      | _ -> assert false
  in
  loop [] xs ys

module GO_term_set = struct
  (* invariant: list of sorted lists *)
  type t = Union of Ontology.PKey.t list list

  let to_sorted_list (Union u) =
    match u with
    | [] -> []
    | [xs] -> xs
    | xs ->
      List.reduce_exn xs ~f:(sorted_list_union ~compare:Ontology.PKey.compare)

  (* let rec sort_lists_by_head = function *)
  (*   | [] -> [] *)
  (*   | [] :: t -> sort_lists_by_head t *)
  (*   | (h1 :: _ as l1) :: t -> *)
  (*     match sort_lists_by_head t with *)
  (*     | [] -> [ l1 ] *)
  (*     | [] :: _ -> assert false *)
  (*     | (h2 :: _ as l2) :: t as sorted_t -> *)
  (*       if Ontology.PKey.compare h1 h2 <= 0 then l1 :: sorted_t *)
  (*       else l2 :: (l1 :: t) *)

  let rec traverse u ~f =
    match find_next u with
    | None -> ()
    | Some (e, u') ->
      f e ; traverse u' ~f
  and find_next = function
    | [] -> None
    | [] :: rest -> find_next rest
    | (h :: t as l0) :: rest ->
      match find_better_candidate h rest with
      | None, rest' -> Some (h, t :: rest')
      | Some h', rest' -> Some (h', l0 :: rest')
  and find_better_candidate x = function
    | [] -> None, []
    | [] :: rest -> find_better_candidate x rest
    | (h :: t as l) :: rest ->
      match Ontology.PKey.compare x h with
      | -1 -> (
          let ans, rest' = find_better_candidate x rest in
          ans, l :: rest'
        )
      | 0 -> (
          let ans, rest' = find_better_candidate x rest in
          match ans with
          | None -> None, t :: rest'
          | Some _ -> ans, l :: rest'
        )
      | 1 -> (
          let ans, rest' = find_better_candidate h rest in
          match ans with
          | None -> Some h, t :: rest'
          | Some _ -> ans, l :: rest'
        )
      | _ -> assert false

  let rec traverse2 xs ys ~f =
    match xs, ys with
    | [], _ -> List.iter ys ~f
    | _, [] -> List.iter xs ~f
    | h_x :: t_x, h_y :: t_y ->
      match Ontology.PKey.compare h_x h_y with
      | -1 -> f h_x ; traverse2 t_x ys ~f
      |  0 -> f h_x ; traverse2 t_x t_y ~f
      |  1 -> f h_y ; traverse2 xs t_y ~f
      | _ -> assert false

  let iter (Union u) ~f =
    match u with
    | [] -> ()
    | [xs] -> List.iter xs ~f
    | [xs ; ys] -> traverse2 xs ys ~f
    | _ -> traverse u ~f
end

let go_categories_by_element ~(element_coordinates:Genomic_interval_collection.t) ~(regulatory_domains:Genomic_interval_collection.t) ~(functional_annot:Functional_annotation.t) =
  (*regulatory domains were constructed for genes that have at least one GO annotation*)
  (*all their IDs should be in the functional annotation*)
  let intersection = Genomic_interval_collection.intersect element_coordinates regulatory_domains in (*dictionary, element ids -> list of gene symbols*)
  List.Assoc.map intersection ~f:(fun neighboors ->
      let cat_lists =
        List.map neighboors ~f:(fun n ->
            let s = Genomic_interval.id n in
            Functional_annotation.extract_terms_exn functional_annot (`Symbol s)
          )
      in
      GO_term_set.Union cat_lists
    )

let elements_by_go_category unique_gocat_by_element =
  let t = String.Table.create () in
  List.iter unique_gocat_by_element ~f:(fun (elt, go_cats) ->
      List.iter go_cats ~f:(fun go_cat -> String.Table.add_multi t ~key:go_cat ~data:(Genomic_interval.id elt))
    ) ;
  String.Table.to_alist t
  |> String.Map.of_alist_exn

let symbol_elements ~(element_coordinates:Genomic_interval_collection.t) ~(regulatory_domains:Genomic_interval_collection.t) =
   let intersection = Genomic_interval_collection.intersect element_coordinates regulatory_domains in (*dictionary, element ids -> list of gene symbols*)
   intersection
