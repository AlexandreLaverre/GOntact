(lang dune 3.1)
(name gontact)
(generate_opam_files true)
(license CECILL-2.1)
(authors "Anamaria Necsulea" "Alexandre Laverré" "Eric Tannier" "Philippe Veber")
(maintainers "Anamaria Necsulea" "Philippe Veber")
(source (uri "https://gitlab.in2p3.fr/anamaria.necsulea/GOntact.git"))
(bug_reports https://gitlab.in2p3.fr/anamaria.necsulea/GOntact/-/issues)
(homepage https://gitlab.in2p3.fr/anamaria.necsulea/GOntact)

(package
  (name gontact)
  (tags (bioinformatics))
  (depends
    (core (>= 0.15.0))
    core_unix))
