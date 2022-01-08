#!/bin/bash

export sp=$1
export dataset=$2
export background=$3
export method=$4

####################################################################################

export path=/beegfs/data/necsulea/GOntact
export pathGO=${path}/data/GeneOntology
export pathEnhancers=${path}/data/enhancers
export pathResults=${path}/results/GO_enrichment_GREAT
export pathScripts=${path}/scripts/GO_enrichment_GREAT

export ensrelease=94

####################################################################################

if [ -e ${pathResults}/${sp}/${method}/${dataset} ]; then
    echo "path output exists"
else
    mkdir -p ${pathResults}/${sp}/${method}/${dataset} 
fi

####################################################################################

perl ${pathScripts}/compute.observed.values.pl --pathInputElements=${pathEnhancers}/${sp}/${dataset}.bed --pathBackgroundElements=${pathEnhancers}/${sp}/${background}.bed --pathGOCategories=${pathGO}/GOCategories.txt --pathGOAnnotations=${pathGO}/${sp}.gene.annotation.txt --pathRegulatoryRegions=${pathResults}/${sp}/${method}/regulatory_regions_Ensembl${ensrelease}.txt --pathOutput=${pathResults}/${sp}/${method}/${dataset}/observed_values_Ensembl${ensrelease}_background${background}.txt 

####################################################################################

## do the test

Rscript ${pathScripts}/test.enrichment.R ${sp} ${dataset} ${background} ${method} > ${pathScripts}/log/test.enrichment.${sp}.${dataset}.${method}.Rout

####################################################################################
