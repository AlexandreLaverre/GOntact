#!/bin/bash

export sp=$1
export genome=$2
export acc=$3

##########################################################

pathResults=../../results/motif_enrichment/${sp}/${acc}
pathGenome=../../data/genome_sequences/${sp}
pathEnhancers=../../data/enhancers/${sp}

##########################################################

for dir in GREAT_only GOntact_only shared
do

    findMotifsGenome.pl ${pathResults}/${dir}/associated_enhancers.bed ${pathGenome}/${genome}.fa ${pathResults}/${dir} -bg ${pathEnhancers}/ENCODE.Laverre2022.bed

done

##########################################################
