#!/bin/bash

###########################################################################################################################

## foreground mouse Vista limb enhancers
## background ENCODE
## GREAT max distance 2Mb

if [ -e results/mouse/Vista_limb_vs_ENCODE/biological_process/GREAT_upstream5kb_downstream1kb_extend2Mb ]; then
    echo "outdir exists"
else
    mkdir -p results/mouse/Vista_limb_vs_ENCODE/biological_process/GREAT_upstream5kb_downstream1kb_extend2Mb 
fi

_build/install/default/bin/gontact --mode=GREAT \
				   --foreground=data/enhancers/human/VistaEnhancers_heart_hg38.bed --background=data/enhancers/human/ENCODE.Laverre2022.bed \
				   --functional-annot=data/GeneOntology/goa_human.gaf  --ontology=data/GeneOntology/go-basic.obo \
				   --gene-annot=data/ensembl_annotations/human/GeneAnnotation_BioMart_Ensembl102_hg38.txt  \
				   --chr-sizes=data/ensembl_annotations/human/chr_sizes_hg38.txt \
				   --upstream=5000 --downstream=1000 --extend=2000000 \
				   --write-foreground --write-background \
				   --output-dir=results/human/Vista_heart_vs_ENCODE/biological_process/GREAT_upstream5kb_downstream1kb_extend2Mb 

## GREAT max distance 1Mb

if [ -e results/human/Vista_heart_vs_ENCODE/biological_process/GREAT_upstream5kb_downstream1kb_extend1Mb ]; then
    echo "outdir exists"
else
    mkdir -p results/human/Vista_heart_vs_ENCODE/biological_process/GREAT_upstream5kb_downstream1kb_extend1Mb 
fi

_build/install/default/bin/gontact --mode=GREAT \
				   --foreground=data/enhancers/human/VistaEnhancers_heart_hg38.bed --background=data/enhancers/human/ENCODE.Laverre2022.bed \
				   --functional-annot=data/GeneOntology/goa_human.gaf  --ontology=data/GeneOntology/go-basic.obo \
				   --gene-annot=data/ensembl_annotations/human/GeneAnnotation_BioMart_Ensembl102_hg38.txt  \
				   --chr-sizes=data/ensembl_annotations/human/chr_sizes_hg38.txt \
				   --upstream=5000 --downstream=1000 --extend=1000000 \
				   --write-foreground --write-background \
				   --output-dir=results/human/Vista_heart_vs_ENCODE/biological_process/GREAT_upstream5kb_downstream1kb_extend1Mb 

## prepare ibed files for contacts

export paths_ibed=""

for file in `ls  data/PCHi-C/human/ibed_files/`
do
    export paths_ibed=data/PCHi-C/human/ibed_files/${file},${paths_ibed}
done

export paths_ibed=${paths_ibed::-1}

## contacts, max distance 1Mb

if [ -e results/human/Vista_heart_vs_ENCODE/biological_process/contacts_mindist25kb_maxdist1Mb ]; then
    echo "outdir exists"
else
    mkdir -p results/human/Vista_heart_vs_ENCODE/biological_process/contacts_mindist25kb_maxdist1Mb 
fi

_build/install/default/bin/gontact --mode=contacts \
				   --foreground=data/enhancers/human/VistaEnhancers_heart_hg38.bed --background=data/enhancers/human/ENCODE.Laverre2022.bed \
				   --functional-annot=data/GeneOntology/goa_human.gaf  --ontology=data/GeneOntology/go-basic.obo \
				   --gene-annot=data/ensembl_annotations/human/GeneAnnotation_BioMart_Ensembl102_hg38.txt  \
				   --chr-sizes=data/ensembl_annotations/human/chr_sizes_hg38.txt \
				   --min-dist-contacts=25000 --max-dist-contacts=1000000 \
				   --ibed-path=${paths_ibed} --bait-coords=data/PCHi-C/human/hg38.baitmap \
				   --write-foreground --write-background \
				   --output-dir=results/human/Vista_heart_vs_ENCODE/biological_process/contacts_mindist25kb_maxdist1Mb


## contacts, max distance 2Mb

if [ -e results/human/Vista_heart_vs_ENCODE/biological_process/contacts_mindist25kb_maxdist2Mb ]; then
    echo "outdir exists"
else
    mkdir -p results/human/Vista_heart_vs_ENCODE/biological_process/contacts_mindist25kb_maxdist2Mb 
fi

_build/install/default/bin/gontact --mode=contacts \
				   --foreground=data/enhancers/human/VistaEnhancers_heart_hg38.bed --background=data/enhancers/human/ENCODE.Laverre2022.bed \
				   --functional-annot=data/GeneOntology/goa_human.gaf  --ontology=data/GeneOntology/go-basic.obo \
				   --gene-annot=data/ensembl_annotations/human/GeneAnnotation_BioMart_Ensembl102_hg38.txt  \
				   --chr-sizes=data/ensembl_annotations/human/chr_sizes_hg38.txt \
				   --min-dist-contacts=25000 --max-dist-contacts=2000000 \
				   --ibed-path=${paths_ibed} --bait-coords=data/PCHi-C/human/hg38.baitmap \
				   --write-foreground --write-background \
				   --output-dir=results/human/Vista_heart_vs_ENCODE/biological_process/contacts_mindist25kb_maxdist2Mb

## hybrid, max distance 1Mb

if [ -e results/human/Vista_heart_vs_ENCODE/biological_process/hybrid_mindist25kb_maxdist1Mb ]; then
    echo "outdir exists"
else
    mkdir -p results/human/Vista_heart_vs_ENCODE/biological_process/hybrid_mindist25kb_maxdist1Mb 
fi

_build/install/default/bin/gontact --mode=hybrid \
				   --foreground=data/enhancers/human/VistaEnhancers_heart_hg38.bed --background=data/enhancers/human/ENCODE.Laverre2022.bed \
				   --functional-annot=data/GeneOntology/goa_human.gaf  --ontology=data/GeneOntology/go-basic.obo \
				   --gene-annot=data/ensembl_annotations/human/GeneAnnotation_BioMart_Ensembl102_hg38.txt  \
				   --chr-sizes=data/ensembl_annotations/human/chr_sizes_hg38.txt \
				   --min-dist-contacts=25000 --max-dist-contacts=1000000 \
				   --upstream=25000 --downstream=25000 --extend=25000 \
				   --ibed-path=${paths_ibed} --bait-coords=data/PCHi-C/human/hg38.baitmap \
				   --write-foreground --write-background \
				   --output-dir=results/human/Vista_heart_vs_ENCODE/biological_process/hybrid_mindist25kb_maxdist1Mb


###########################################################################################################################

_build/install/default/bin/gontact --mode=contacts \
				   --foreground=data/enhancers/mouse/VistaLimbEnhancerGenie_mm10.bed --background=data/enhancers/mouse/ENCODE.Laverre2022.bed \
				   --functional-annot=data/GeneOntology/mgi.gaf  --ontology=data/GeneOntology/go-basic.obo \
				   --gene-annot=data/ensembl_annotations/mouse/GeneAnnotation_BioMart_Ensembl102_mm10.txt \
				   --chr-sizes=data/ensembl_annotations/mouse/chr_sizes_mm10.txt \
				   --min-dist-contacts=25000 --max-dist-contacts=1000000 \
				   --ibed-path=${paths_ibed} --bait-coords=data/PCHi-C/mouse/mm10.baitmap \
				   --write-foreground --write-background \
				   --output-dir= results/mouse/VistaLimbEnhancerGenie_vs_ENCODE/biological_processcontacts_mindist25kb_maxdist1Mb

