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

## GREAT max distance 2Mb

if [ -e results/${sp}/${outprefix}/biological_process/GREAT_upstream5kb_downstream1kb_extend2Mb ]; then
    echo "outdir exists"
else
    mkdir -p results/${sp}/${outprefix}/biological_process/GREAT_upstream5kb_downstream1kb_extend2Mb 
fi

_build/install/default/bin/gontact --mode=GREAT \
				   --foreground=data/enhancers/${sp}/${foreground}.bed --background=data/enhancers/${sp}/${background}.bed \
				   --functional-annot=data/GeneOntology/${goa_file}.gaf  --ontology=data/GeneOntology/go-basic.obo \
				   --gene-annot=data/ensembl_annotations/${sp}/GeneAnnotation_BioMart_Ensembl102_${genome}.txt  \
				   --chr-sizes=data/ensembl_annotations/${sp}/chr_sizes_${genome}.txt \
				   --upstream=5000 --downstream=1000 --extend=2000000 \
				   --write-foreground --write-background \
				   --output-dir=results/${sp}/${outprefix}/biological_process/GREAT_upstream5kb_downstream1kb_extend2Mb 

## GREAT max distance 1Mb

if [ -e results/${sp}/${outprefix}/biological_process/GREAT_upstream5kb_downstream1kb_extend1Mb ]; then
    echo "outdir exists"
else
    mkdir -p results/${sp}/${outprefix}/biological_process/GREAT_upstream5kb_downstream1kb_extend1Mb 
fi

_build/install/default/bin/gontact --mode=GREAT \
				   --foreground=data/enhancers/${sp}/${foreground}.bed --background=data/enhancers/${sp}/${background}.bed \
				   --functional-annot=data/GeneOntology/${goa_file}.gaf  --ontology=data/GeneOntology/go-basic.obo \
				   --gene-annot=data/ensembl_annotations/${sp}/GeneAnnotation_BioMart_Ensembl102_${genome}.txt  \
				   --chr-sizes=data/ensembl_annotations/${sp}/chr_sizes_${genome}.txt \
				   --upstream=5000 --downstream=1000 --extend=1000000 \
				   --write-foreground --write-background \
				   --output-dir=results/${sp}/${outprefix}/biological_process/GREAT_upstream5kb_downstream1kb_extend1Mb 

###########################################################################################################################

## prepare ibed files for contacts

export paths_ibed=""

for file in `ls  data/PCHi-C/${sp}/ibed_files/`
do
    export paths_ibed=data/PCHi-C/${sp}/ibed_files/${file},${paths_ibed}
done

## chop off the comma
export paths_ibed=${paths_ibed::-1}

###########################################################################################################################

## contacts, max distance 1Mb

if [ -e results/${sp}/${outprefix}/biological_process/contacts_mindist25kb_maxdist1Mb ]; then
    echo "outdir exists"
else
    mkdir -p results/${sp}/${outprefix}/biological_process/contacts_mindist25kb_maxdist1Mb 
fi

_build/install/default/bin/gontact --mode=contacts \
				   --foreground=data/enhancers/${sp}/${foreground}.bed --background=data/enhancers/${sp}/${background}.bed \
				   --functional-annot=data/GeneOntology/${goa_file}.gaf  --ontology=data/GeneOntology/go-basic.obo \
				   --gene-annot=data/ensembl_annotations/${sp}/GeneAnnotation_BioMart_Ensembl102_${genome}.txt  \
				   --chr-sizes=data/ensembl_annotations/${sp}/chr_sizes_${genome}.txt \
				   --min-dist-contacts=25000 --max-dist-contacts=1000000 \
				   --ibed-path=${paths_ibed} --bait-coords=data/PCHi-C/${sp}/${genome}.baitmap \
				   --write-foreground --write-background \
				   --output-dir=results/${sp}/${outprefix}/biological_process/contacts_mindist25kb_maxdist1Mb


## contacts, max distance 2Mb

if [ -e results/${sp}/${outprefix}/biological_process/contacts_mindist25kb_maxdist2Mb ]; then
    echo "outdir exists"
else
    mkdir -p results/${sp}/${outprefix}/biological_process/contacts_mindist25kb_maxdist2Mb 
fi

_build/install/default/bin/gontact --mode=contacts \
				   --foreground=data/enhancers/${sp}/${foreground}.bed --background=data/enhancers/${sp}/${background}.bed \
				   --functional-annot=data/GeneOntology/${goa_file}.gaf  --ontology=data/GeneOntology/go-basic.obo \
				   --gene-annot=data/ensembl_annotations/${sp}/GeneAnnotation_BioMart_Ensembl102_${genome}.txt  \
				   --chr-sizes=data/ensembl_annotations/${sp}/chr_sizes_${genome}.txt \
				   --min-dist-contacts=25000 --max-dist-contacts=2000000 \
				   --ibed-path=${paths_ibed} --bait-coords=data/PCHi-C/${sp}/${genome}.baitmap \
				   --write-foreground --write-background \
				   --output-dir=results/${sp}/${outprefix}/biological_process/contacts_mindist25kb_maxdist2Mb

###########################################################################################################################

## hybrid, max distance 1Mb

if [ -e results/${sp}/${outprefix}/biological_process/hybrid_mindist25kb_maxdist1Mb ]; then
    echo "outdir exists"
else
    mkdir -p results/${sp}/${outprefix}/biological_process/hybrid_mindist25kb_maxdist1Mb 
fi

_build/install/default/bin/gontact --mode=hybrid \
				   --foreground=data/enhancers/${sp}/${foreground}.bed --background=data/enhancers/${sp}/${background}.bed \
				   --functional-annot=data/GeneOntology/${goa_file}.gaf  --ontology=data/GeneOntology/go-basic.obo \
				   --gene-annot=data/ensembl_annotations/${sp}/GeneAnnotation_BioMart_Ensembl102_${genome}.txt  \
				   --chr-sizes=data/ensembl_annotations/${sp}/chr_sizes_${genome}.txt \
				   --min-dist-contacts=25000 --max-dist-contacts=1000000 \
				   --upstream=25000 --downstream=25000 --extend=25000 \
				   --ibed-path=${paths_ibed} --bait-coords=data/PCHi-C/${sp}/${genome}.baitmap \
				   --write-foreground --write-background \
				   --output-dir=results/${sp}/${outprefix}/biological_process/hybrid_mindist25kb_maxdist1Mb

###########################################################################################################################
