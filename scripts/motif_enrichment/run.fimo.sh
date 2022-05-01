#!/bin/bash

export sp=$1
export motif=$2
export acc=$3
export dir=$4


##########################################################

pathJASPAR=../../data/JASPAR
pathResults=../../results/motif_enrichment/${sp}

##########################################################

fimo -o ${pathResults}/${acc}/${dir}/${motif}  ${pathJASPAR}/${motif}.meme ${pathResults}/${acc}/${dir}/associated_enhancers.fa

##########################################################
