#!/bin/bash

export sp=$1
export dataset=$2
export method=$3

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

perl ${pathScripts}/compute.observed.values.pl --pathInputElements=${pathEnhancers}/${sp}/${dataset}.bed --pathGOCategories=${pathGO}/GOCategories.txt --pathGOAnnotations=${pathGO}/${sp}.gene.annotation.txt --pathRegulatoryRegions=${pathResults}/${sp}/${method}/regulatory_regions_Ensembl${ensrelease}.txt --pathOutput=${pathResults}/${sp}/${method}/${dataset}/observed_values_Ensembl${ensrelease}.txt 

####################################################################################

## do the test

R CMD BATCH '--args sp=\"${sp}\" dataset=\"${dataset}\" method=\"${method}\"' ${pathScripts}/test.enrichment.R log/${pathScripts}/test.enrichment.${sp}.${dataset}.${method}.Rout

####################################################################################
