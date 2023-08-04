open Core

type t = {
  bait_chr : string ;
  bait_start : int ;
  bait_end : int ;
  bait_name : string ;
  otherEnd_chr : string ;
  otherEnd_start : int ;
  otherEnd_end : int ;
  otherEnd_name : string ;
  n_reads : int ;
  score : float ;
}
[@@deriving sexp]


let distance cc =
  match String.equal cc.bait_chr cc.otherEnd_chr with
  | true -> (
      let midpos1 = (cc.bait_start + cc.bait_end) / 2 in
      let midpos2 = (cc.otherEnd_start + cc.otherEnd_end) / 2 in
      let dist = Int.abs (midpos1 - midpos2) in
      Some dist
    )
  | false -> None

let contacted_fragment i =
  let (chr, start_pos, end_pos) = (i.otherEnd_chr, i.otherEnd_start, i.otherEnd_end) in
  Genomic_interval.make chr start_pos end_pos Unstranded

let get_id_frag i =
  let (chr, start_pos, end_pos) = (i.otherEnd_chr, i.otherEnd_start, i.otherEnd_end) in
  Printf.sprintf "%s:%d-%d" chr start_pos end_pos

let get_id_bait i =
  let (chr, start_pos, end_pos) = (i.bait_chr, i.bait_start, i.bait_end) in
  Printf.sprintf "%s:%d-%d" chr start_pos end_pos

let compare i1 i2 =
  match String.compare i1.bait_chr i2.bait_chr with
  | 0 -> (
      match Int.compare i1.bait_start i2.bait_start with
      | 0 -> (
          match Int.compare i1.bait_end i2.bait_end with
          | 0 -> (
              match String.compare i1.otherEnd_chr i2.otherEnd_chr with
              | 0 -> (
                  match Int.compare i1.otherEnd_start i2.otherEnd_start with
                  | 0 -> Int.compare i1.otherEnd_end i2.otherEnd_end
                  | n -> n
                )
              | n -> n
            )
          |  n -> n
        )
      | n -> n
    )
  | n -> n

(*Stdlib.compare (i1.bait_chr, i1.bait_start, i1.bait_end, i1.otherEnd_chr, i1.otherEnd_start, i1.otherEnd_end) (i2.bait_chr, i2.bait_start, i2.bait_end, i2.otherEnd_chr, i2.otherEnd_start, i2.otherEnd_end) *)

(*  let id1 = Printf.sprintf "%s %s" (get_id_bait i1) (get_id_frag i1) in
  let id2 = Printf.sprintf "%s %s" (get_id_bait i2) (get_id_frag i2) in
  String.compare id1 id2
*)

let equal i1 i2 = (compare i1 i2 = 0)
  (*let id1 = Printf.sprintf "%s %s" (get_id_bait i1) (get_id_frag i1) in
  let id2 = Printf.sprintf "%s %s" (get_id_bait i2) (get_id_frag i2) in
  String.equal id1 id2
  *)

let hash x = Hashtbl.hash (x.bait_chr,x.bait_start, x.bait_end, x.otherEnd_chr, x.otherEnd_start, x.otherEnd_end)

let symbol_annotate_baits ~bait_collection ~genome_annotation ~max_dist =
  let tss_intervals = Genomic_annotation.all_tss_intervals genome_annotation max_dist in
  let intersection = (*String.Map - key = bait ids ; values = list of gene_id:gene_symbol mixed id*)
    Genomic_interval_collection.intersect bait_collection tss_intervals
    |> List.map ~f:(fun (elt, fragments) -> Genomic_interval.id elt, List.map fragments ~f:Genomic_interval.id)
    |> String.Map.of_alist_exn
  in
  let get_symbol id =
    match String.split id ~on:':' with
    | [ _ ; symbol ] -> Some symbol
    | _ -> None (*not optimal, fix this*)
  in
  let symbol_annot = String.Map.map intersection ~f:(fun l -> List.filter_map l ~f:get_symbol) in
  String.Map.map symbol_annot ~f:(fun l -> List.dedup_and_sort ~compare:String.compare l)

let output_bait_annotation ~bait_collection ~bait_annotation ~path =
  let bait_list = Genomic_interval_collection.interval_list bait_collection in
  let id_baits = List.map bait_list ~f:(fun b -> Genomic_interval.id b) in
  Out_channel.with_file path ~append:false ~f:(fun output ->
      Out_channel.output_string output "BaitID\tGOID\n" ;
      List.iter id_baits  ~f:(fun id ->
          match (String.Map.find bait_annotation id) with
          | None ->  Printf.fprintf output "%s\t\n" id
          | Some l -> Printf.fprintf output "%s\t%s\n" id (String.concat ~sep:"," l)
        )
    )

let annotations_by_element ~(element_coordinates:Genomic_interval_collection.t) ~(fragments:Genomic_interval_collection.t) ~fragment_to_baits ~annotated_baits =
  (* for each element, find the list of annotations - gene symbols or GO categories - associated with it *)
  let intersection =
    Genomic_interval_collection.intersect element_coordinates fragments
    |> List.map ~f:(fun (elt, fragments) -> Genomic_interval.id elt, List.map fragments ~f:Genomic_interval.id)
    |> String.Map.of_alist_exn
  in
  let elbaits =
    String.Map.map intersection ~f:(fun l ->
        List.filter_map l ~f:(fun frag -> String.Map.find fragment_to_baits frag)
        |> List.join
        |> List.dedup_and_sort ~compare:String.compare
      )
  in
  String.Map.map elbaits ~f:(fun l ->
      List.filter_map l ~f:(fun bait -> String.Map.find annotated_baits bait)
      |> List.join
      |> List.dedup_and_sort ~compare:String.compare
    )

let elements_by_annotation elgolist =
  let tuples =
    List.map elgolist ~f:(fun (elt, l) ->
        List.map l ~f:(fun g -> (g, Genomic_interval.id elt))
      )
    |> List.join
  in
  let gomap = String.Map.of_alist_multi tuples in
  gomap

let write_contact contact output =
  Printf.fprintf output "%s\t%d\t%d\t%s\t%s\t%d\t%d\t%s\t%d\t%f\n" contact.bait_chr contact.bait_start contact.bait_end contact.bait_name  contact.otherEnd_chr contact.otherEnd_start contact.otherEnd_end contact.otherEnd_name 0 0.0
