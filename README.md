# GOntact

**GOntact** is a tool for functional enrichment analysis of *cis*-regulatory elements using chromatin contact data. Unlike traditional proximity-based approaches, GOntact leverages promoter-capture Hi-C (PCHi-C) data to infer gene-enhancer relationships and derive Gene Ontology enrichments. 

It implements the method described in [Laverré et al., 2022 (bioRxiv)](https://www.biorxiv.org/content/10.1101/2022.06.13.495495v1), and uses PCHi-C data from [Laverré et al., Genome Research, 2022](https://pmc.ncbi.nlm.nih.gov/articles/PMC8805723/) to provide biologically coherent functional annotations and novel insights into gene regulation across human and mouse genomes.

[GOntact installation](#install)

[GOntact usage](#usage)

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
opam install dune
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
that follow are for the `gontact` executable rather than for `dune): 

```
dune exec gontact -- --help
```

All GOntact commands thereafter will start with `dune exec gontact --
` instead of `gontact`. 

### Installing the GOntact package directly using opam 


`GOntact` can be installed via [opam](http://opam.ocaml.org/), using
the following command:

```
opam pin add -y gontact https://gitlab.in2p3.fr/anamaria.necsulea/GOntact.git
```

This will install an executable named `gontact` in your .opam directory. This executable
should be accessible from any directory, if your PATH variable
was correctly configured when installing opam.

## Usage
<a name="usage"></a>

For the following usage commands, the GOntact executable is named
`gontact`. Please refer to the above instructions if you
installed GOntact directly from source.

To test GOntact, you can use this example dataset (make sure you don't
have a `data` directory where you run this!):

```
wget http://pbil.univ-lyon1.fr/members/necsulea/GOntact/data.tar.gz
test ! -d data && tar -xzvf data.tar.gz
```

These commands will create a directory named `data`, containing
several data types.

### Data types
[Chromatin contact data](#contacts)

[Gene annotations](#annot)

[Gene Ontology annotations](#GO)

[Enhancer coordinates](#enhancers)

#### Chromatin contact data
<a name="contacts"></a>
This dataset includes chromatin contacts detected with the
PCHi-C approach, for human and mouse. This data was described in
[Laverré et al., Genome Research, 2022](https://pmc.ncbi.nlm.nih.gov/articles/PMC8805723/). Chromatin interactions were
scored for several cell types for each species and are provided in the
[ibed format](https://www.bioconductor.org/packages/devel/bioc/vignettes/Chicago/inst/doc/Chicago.html),
in the `data/PCHi-C subfolder`. 

Note that these chromatin contacts are provided with respect to the
GRCh38 (hg38) genome assembly for human and to the GRCm38 (mm10)
genome assembly for mouse. 

#### Gene annotations
<a name="annot"></a>
GOntact also needs a set of genomic annotations as an
input. These should be provided in the [GTF format](https://www.ensembl.org/info/website/upload/gff.html).  Example files for human and mouse are provided in the subfolder `data/ensembl_annotations`. They were downloaded from the Ensembl database and correspond to the GRCh38 (hg38) genome assembly for human and to the GRCm38 (mm10)
genome assembly for mouse. 

#### Gene Ontology annotation
<a name="GO"></a>
Gene Ontology annotations are provided for human and mouse in the
`data/GeneOntology` subfolder (goa_human.gaf for human and mgi.gaf for
mouse). These files correspond to the 2025-07-22 release of the Gene Ontology database. They were downloaded from
geneontology.org. 

#### Enhancer coordinates
<a name="enhancers"></a>
Example enhancer datasets are also provided in the subfolder
`data/enhancers`.


### Basic usage
Here is an example of a command line that runs GOntact in "GREAT" mode,
using Vista midbrain enhancers as a foreground set and the full set of
ENCODE enhancers as a background set:

```
mkdir GREAT_results_bp

gontact enrich \
--mode=GREAT \
--gene-annot=data/genomic_annotations/Homo_sapiens.GRCh38.115.gtf \
--functional-annot=data/GeneOntology/goa_human.gaf \
--ontology=data/GeneOntology/go-basic.obo \
--foreground=data/enhancers/human/VistaEnhancers_midbrain_hg38.bed \
--background=data/enhancers/human/ENCODE.Laverre2022.bed \
--upstream=5000 \
--downstream=1000 \
--extend=1000000 \
--output-dir=GREAT_results_bp \
--output-prefix=VistaEnhancers_midbrain \
--domain=biological_process 
```

Here is an example of a command line that runs GOntact in "contacts" mode,
using Vista midbrain enhancers as a foreground set, the full set of
ENCODE enhancers as a background set, and the set of PCHi-C
interactions that are shared in at least two samples.

```
mkdir contacts_results_bp

gontact enrich \
--mode=contacts \
--gene-annot=data/ensembl_annotations/Homo_sapiens.GRCh38.115.gtf \
--functional-annot=data/GeneOntology/goa_human.gaf \
--ontology=data/GeneOntology/go-basic.obo \
--foreground=data/enhancers/human/VistaEnhancers_midbrain_hg38.bed \
--background=data/enhancers/human/ENCODE.Laverre2022.bed \
--ibed-path=data/PCHi-C/human/ibed_files/shared_contacts_min2samples.ibed \
--min-score=0 \
--bait-coords=data/PCHi-C/human/hg38.baitmap \
--output-dir=contacts_results_bp \
--output-prefix=VistaEnhancers_midbrain \
--domain=biological_process
```
