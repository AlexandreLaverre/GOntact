#!/bin/bash

export sp=$1
export genome=$2
export acc=$3
export dir=$4

##########################################################

pathResults=../../results/motif_enrichment/${sp}
pathGenome=../../data/genome_sequences/${sp}

##########################################################

findMotifsGenome.pl ${pathResults}/${acc}/${dir}/associated_enhancers.bed ${pathGenome}/${genome}.fa ${pathResults}/${acc}/${dir} -bg ${pathResults}/ENCODE.Laverre2022.random50K.bed

##########################################################
