#!/bin/bash

export sp=$1

####################################################################################

export path=/beegfs/data/necsulea/GOntact
export pathGO=${path}/data/GeneOntology
export pathResults=${path}/results/GO_enrichment_GREAT/${sp}
export pathScripts=${path}/scripts/GO_enrichment_GREAT

export ensrelease=94
export basal5=5000
export basal3=1000
export maxextend=1e+06

####################################################################################

perl ${pathScripts}/compute.expected.values.pl --pathGOCategories=${pathGO}/GOCategories.txt --pathGOAnnotations=${pathGO}/${sp}.gene.annotation.txt --pathRegulatoryRegions=${pathResults}/${sp}/regulatory_regions_Ensembl${ensrelease}_basal5${basal5}_basal3${basal3}_max_extend${maxextend}.txt --pathOutput=${pathResults}/${sp}/expected_values_Ensembl${ensrelease}_basal5${basal5}_basal3${basal3}_max_extend${maxextend}.txt 

####################################################################################
