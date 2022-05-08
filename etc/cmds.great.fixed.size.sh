#!/bin/bash

export sp=$1
export foreground=$2
export background=$3
export outprefix=$4

############################################################################################

if [ ${sp} = "human" ]; then
    export goa_file=goa_human
    export genome=hg38
fi

if [ ${sp} = "mouse" ]; then
    export goa_file=mgi
    export genome=mm10
fi

############################################################################################

## GREAT fixed size max distance 1Mb

if [ -e results/${sp}/${outprefix}/biological_process/GREAT_fixed_size_upstream1Mb_downstream1Mb ]; then
    echo "outdir exists"
else
    mkdir -p results/${sp}/${outprefix}/biological_process/GREAT_fixed_size_upstream1Mb_downstream1Mb 
fi

_build/install/default/bin/gontact --mode=GREAT \
				   --foreground=data/enhancers/${sp}/${foreground}.bed --background=data/enhancers/${sp}/${background}.bed \
				   --functional-annot=data/GeneOntology/${goa_file}.gaf  --ontology=data/GeneOntology/go-basic.obo \
				   --gene-annot=data/ensembl_annotations/${sp}/GeneAnnotation_BioMart_Ensembl102_${genome}.txt  \
				   --chr-sizes=data/ensembl_annotations/${sp}/chr_sizes_${genome}.txt \
				   --upstream=1000000 --downstream=1000000 --extend=0 \
				   --write-foreground --write-background \
				   --output-dir=results/${sp}/${outprefix}/biological_process/GREAT_fixed_size_upstream1Mb_downstream1Mb

###########################################################################################################################
