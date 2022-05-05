open Core

let parse_line l =
  match String.split l ~on:'\t' with
  | [ go_id; enhancer_ids ] ->
    let enhancer_list = String.split enhancer_ids ~on:',' in
    List.map enhancer_list ~f:(fun enh_id -> (enh_id, go_id))
  | _ -> invalid_arg "invalid enhancer annotation file"

let from_gontact_result_file path =
 let lines = In_channel.read_lines path in 
 let all_annot = List.concat_map lines ~f:parse_line in
 String.Map.of_alist_multi all_annot

let common_annotations annot1 annot2 =
  let keys1 = String.Set.of_list (String.Map.keys annot1) in
  let keys2 = String.Set.of_list (String.Map.keys annot2) in
  let common_keys = String.Set.to_list (String.Set.inter keys1 keys2) in
  let nb_common_list = List.map common_keys ~f:(fun id -> (
        let go1 = String.Map.find_exn annot1 id in
        let go2 = String.Map.find_exn annot2 id in
        let nb_common = String.Set.length (String.Set.inter (String.Set.of_list go1) (String.Set.of_list go2)) in
        (id, nb_common))) in
  String.Map.of_alist_exn nb_common_list

let compare_annotations annot1 annot2 path =
  let keys1 = String.Set.of_list (String.Map.keys annot1) in
  let keys2 = String.Set.of_list (String.Map.keys annot2) in
  let all_keys = String.Set.to_list (String.Set.union keys1 keys2) in 
  let ca = common_annotations annot1 annot2 in
  Out_channel.with_file path ~append:false ~f:(fun output ->
      Printf.fprintf output "EnhancerID\tNb1\tNb2\tNbComon\n" ;
      List.iter all_keys ~f:(fun key -> (
          let nb1 = (if (String.Map.mem annot1 key) then (List.length (String.Map.find_exn annot1 key)) else 0) in
          let nb2 = (if (String.Map.mem annot2 key) then (List.length (String.Map.find_exn annot2 key)) else 0) in
          let nbc = (if (String.Map.mem ca key) then (String.Map.find_exn ca key) else 0) in
          Printf.fprintf output "%s\t%d\t%d\t%d\n" key nb1 nb2 nbc 
        )))

  
