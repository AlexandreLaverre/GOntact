# GOntact
[Installation](#install)
[Usage](#usage)

## Installation
<a name="install"></a>

### Installing GOntact from source

To install GOntact from source, first clone the GitLab repository
using the following command: 


```
git clone git@gitlab.in2p3.fr:anamaria.necsulea/GOntact.git
```

To compile GOntact, you will need the [opam](http://opam.ocaml.org/)
package manager for OCaml. On Debian-based Linux distributions, you can install it and
initialize it using the following
commands: 

```
sudo apt install opam
opam init
```

You will also need the [dune](https://github.com/ocaml/dune) build
system for OCaml.

```
sudo apt install ocaml-dune
```

Once opam and dune are installed, in the GOntact directory created before, install the required
dependencies using the following commands:

```
opam pin add -n gontact .
opam install gontact --deps-only

opam pin add -n gontact-server .
opam install gontact-server --deps-only
```

You can now compile GOntact with the following command, run from the
GOntact directory: 

```
dune build
```

To use this version of GOntact, you can use the following commands
(note the extra "--", which signifies that the command-line arguments
that follow are for the gontact executable rather than for dune): 

```
dune exec gontact -- --help
```

All GOntact commands thereafter will start with 'dune exec gontact --
' instead of 'gontact'. 

### Installing the GOntact package directly using opam 


`GOntact` can be installed via [opam](http://opam.ocaml.org/), using
the following command:

```
opam pin add -y gontact https://gitlab.in2p3.fr/anamaria.necsulea/GOntact.git
```

## Usage
<a name="usage"></a>

To test GOntact, you can use this example dataset (make sure you don't
have a `data` directory where you run this!):

```
wget http://pbil.univ-lyon1.fr/members/necsulea/GOntact/data.tar.gz
test ! -d data && tar -xzvf data.tar.gz
```

This dataset includes chromatin interactions detected with the
PCHi-C approach, for human and mouse. This data was described in
Laverré et al., Genome Research, 2022. Chromatin interactions were
scored for several cell types for each species and are provided in the
[ibed format](https://www.bioconductor.org/packages/devel/bioc/vignettes/Chicago/inst/doc/Chicago.html),
in the data/PCHi-C subfolder.

GOntact also needs a set of gene annotations as an
input. Annotations were downloaded from the Ensembl database (release
102) and are
provided in the subfolder data/ensembl_annotations. This subfolder
also contains files with chromosome size information, which are needed
to run GOntact in "GREAT" mode.

Gene Ontology annotations are provided for human and mouse in the
data/GeneOntology subfolder (goa_human.gaf for human and mgi.gaf for
mouse). These files were downloaded on 07/12/2021 from
geneontology.org. The Gene Ontology hierarchy

Example enhancer datasets are also provided in the subfolder
data/enhancers.

Here is an example of a command line that runs GOntact in "GREAT" mode,
using Vista midbrain enhancers as a foreground set and the full set of
ENCODE enhancers as a background set:

```
mkdir GREAT_results_bp

dune exec gontact -- \
--mode=GREAT \
--gene-annot=data/ensembl_annotations/human/GeneAnnotation_BioMart_Ensembl102_hg38.txt \
--functional-annot=data/GeneOntology/goa_human.gaf \
--ontology=data/GeneOntology/go-basic.obo \
--foreground=data/enhancers/human/VistaEnhancers_midbrain_hg38.bed \
--background=data/enhancers/human/ENCODE.Laverre2022.bed \
--chr-sizes=data/ensembl_annotations/human/chr_sizes_hg38.txt \
--upstream=5000 \
--downstream=1000 \
--extend=1000000 \
--output-dir=GREAT_results_bp \
--output-prefix=VistaEnhancers_midbrain \
--domain=biological_process \
```


Here is an example of a command line that runs GOntact in "contacts" mode,
using Vista midbrain enhancers as a foreground set, the full set of
ENCODE enhancers as a background set, and the set of PCHi-C
interactions that are shared in at least two samples.

```
mkdir contacts_results_bp

dune exec gontact -- \
--mode=contacts \
--gene-annot=data/ensembl_annotations/human/GeneAnnotation_BioMart_Ensembl102_hg38.txt \
--functional-annot=data/GeneOntology/goa_human.gaf \
--ontology=data/GeneOntology/go-basic.obo \
--foreground=data/enhancers/human/VistaEnhancers_midbrain_hg38.bed \
--background=data/enhancers/human/ENCODE.Laverre2022.bed \
--chr-sizes=data/ensembl_annotations/human/chr_sizes_hg38.txt \
--ibed-path=data/PCHi-C/human/ibed_files/shared_contacts_min2samples.ibed \
--min-score=0 \
--bait-coords=data/PCHi-C/human/hg38.baitmap \
--output-dir=contacts_results_bp \
--output-prefix=VistaEnhancers_midbrain \
--domain=biological_process
```
