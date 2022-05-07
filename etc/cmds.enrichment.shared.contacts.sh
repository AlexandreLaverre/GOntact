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

###########################################################################################################################

export paths_ibed=data/PCHi-C/${sp}/ibed_files/shared_contacts_min2samples.ibed

###########################################################################################################################

## contacts, max distance 1Mb, min distance 25kb

if [ -e results/${sp}/${outprefix}/biological_process/shared_contacts_mindist25kb_maxdist1Mb ]; then
    echo "outdir exists"
else
    mkdir -p results/${sp}/${outprefix}/biological_process/shared_contacts_mindist25kb_maxdist1Mb
fi

_build/install/default/bin/gontact --mode=shared_contacts \
				   --foreground=data/enhancers/${sp}/${foreground}.bed --background=data/enhancers/${sp}/${background}.bed \
				   --functional-annot=data/GeneOntology/${goa_file}.gaf  --ontology=data/GeneOntology/go-basic.obo \
				   --gene-annot=data/ensembl_annotations/${sp}/GeneAnnotation_BioMart_Ensembl102_${genome}.txt  \
				   --chr-sizes=data/ensembl_annotations/${sp}/chr_sizes_${genome}.txt \
				   --min-dist-contacts=25000 --max-dist-contacts=1000000 \
				   --min-score=0 \
				   --ibed-path=${paths_ibed} --bait-coords=data/PCHi-C/${sp}/${genome}.baitmap \
				   --write-foreground --write-background \
				   --output-dir=results/${sp}/${outprefix}/biological_process/shared_contacts_mindist25kb_maxdist1Mb


## contacts, max distance 2Mb, min distance 25kb

if [ -e results/${sp}/${outprefix}/biological_process/shared_contacts_mindist25kb_maxdist2Mb ]; then
    echo "outdir exists"
else
    mkdir -p results/${sp}/${outprefix}/biological_process/shared_contacts_mindist25kb_maxdist2Mb
fi

_build/install/default/bin/gontact --mode=contacts \
				   --foreground=data/enhancers/${sp}/${foreground}.bed --background=data/enhancers/${sp}/${background}.bed \
				   --functional-annot=data/GeneOntology/${goa_file}.gaf  --ontology=data/GeneOntology/go-basic.obo \
				   --gene-annot=data/ensembl_annotations/${sp}/GeneAnnotation_BioMart_Ensembl102_${genome}.txt  \
				   --chr-sizes=data/ensembl_annotations/${sp}/chr_sizes_${genome}.txt \
				   --min-dist-contacts=25000 --max-dist-contacts=2000000 \
				   --min-score=0 \
				   --ibed-path=${paths_ibed} --bait-coords=data/PCHi-C/${sp}/${genome}.baitmap \
				   --write-foreground --write-background \
				   --output-dir=results/${sp}/${outprefix}/biological_process/shared_contacts_mindist25kb_maxdist2Mb

###########################################################################################################################

## contacts, max distance 1Mb, min distance 0 kb

if [ -e results/${sp}/${outprefix}/biological_process/shared_contacts_mindist0kb_maxdist1Mb ]; then
    echo "outdir exists"
else
    mkdir -p results/${sp}/${outprefix}/biological_process/shared_contacts_mindist0kb_maxdist1Mb
fi

_build/install/default/bin/gontact --mode=contacts \
				   --foreground=data/enhancers/${sp}/${foreground}.bed --background=data/enhancers/${sp}/${background}.bed \
				   --functional-annot=data/GeneOntology/${goa_file}.gaf  --ontology=data/GeneOntology/go-basic.obo \
				   --gene-annot=data/ensembl_annotations/${sp}/GeneAnnotation_BioMart_Ensembl102_${genome}.txt  \
				   --chr-sizes=data/ensembl_annotations/${sp}/chr_sizes_${genome}.txt \
				   --min-dist-contacts=0 --max-dist-contacts=1000000 \
				   --min-score=0 \
				   --ibed-path=${paths_ibed} --bait-coords=data/PCHi-C/${sp}/${genome}.baitmap \
				   --write-foreground --write-background \
				   --output-dir=results/${sp}/${outprefix}/biological_process/shared_contacts_mindist0kb_maxdist1Mb


## contacts, max distance 2Mb

if [ -e results/${sp}/${outprefix}/biological_process/shared_contacts_mindist0kb_maxdist2Mb ]; then
    echo "outdir exists"
else
    mkdir -p results/${sp}/${outprefix}/biological_process/shared_contacts_mindist0kb_maxdist2Mb
fi

_build/install/default/bin/gontact --mode=contacts \
				   --foreground=data/enhancers/${sp}/${foreground}.bed --background=data/enhancers/${sp}/${background}.bed \
				   --functional-annot=data/GeneOntology/${goa_file}.gaf  --ontology=data/GeneOntology/go-basic.obo \
				   --gene-annot=data/ensembl_annotations/${sp}/GeneAnnotation_BioMart_Ensembl102_${genome}.txt  \
				   --chr-sizes=data/ensembl_annotations/${sp}/chr_sizes_${genome}.txt \
				   --min-dist-contacts=0 --max-dist-contacts=2000000 \
				   --min-score=0 \
				   --ibed-path=${paths_ibed} --bait-coords=data/PCHi-C/${sp}/${genome}.baitmap \
				   --write-foreground --write-background \
				   --output-dir=results/${sp}/${outprefix}/biological_process/shared_contacts_mindist0kb_maxdist2Mb

###########################################################################################################################
