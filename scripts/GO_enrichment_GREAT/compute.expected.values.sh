#!/bin/bash

export sp=$1
export method=$2 ## classical_upstream5kb_downstream1kb_extend1Mb

####################################################################################

export path=/beegfs/data/necsulea/GOntact
export pathGO=${path}/data/GeneOntology
export pathGenomes=${path}/data/genome_sequences
export pathResults=${path}/results/GO_enrichment_GREAT
export pathScripts=${path}/scripts/GO_enrichment_GREAT

export ensrelease=94

####################################################################################

perl ${pathScripts}/compute.expected.values.pl --pathGOCategories=${pathGO}/GOCategories.txt --pathGOAnnotations=${pathGO}/${sp}.gene.annotation.txt --pathRegulatoryRegions=${pathResults}/${sp}/${method}/regulatory_regions_Ensembl${ensrelease}.txt --pathGenomeSequence=${pathGenomes}/${sp}/genome_Ensembl${ensrelease}.fa --pathOutput=${pathResults}/${sp}/${method}/expected_values_Ensembl${ensrelease}.txt 

####################################################################################
