#!/bin/bash

export sp=$1
export genome=$2
export acc=$3

##########################################################

pathResults=../../results/motif_enrichment/${sp}
pathGenome=../../data/genome_sequences/${sp}

##########################################################

for dir in GREAT_only GOntact_only shared
do
    findMotifsGenome.pl ${pathResults}/${acc}/${dir}/associated_enhancers.bed ${pathGenome}/${genome}.fa ${pathResults}/${acc}/${dir} -bg ${pathResults}/ENCODE.Laverre2022.random50K.bed
done

##########################################################
