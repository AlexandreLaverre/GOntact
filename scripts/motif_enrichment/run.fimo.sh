#!/bin/bash

export sp=$1
export motif=$2
export acc=$3
export dir=$4


##########################################################

pathJASPAR=../../data/JASPAR
pathResults=../../results/motif_enrichment/${sp}

##########################################################

if [ -e ${pathResults}/${acc}/${dir}/${motif} ]; then
    echo "already done"
else
    fimo -o ${pathResults}/${acc}/${dir}/${motif}  ${pathJASPAR}/${motif}.meme ${pathResults}/${acc}/${dir}/associated_enhancers.fa
fi

##########################################################

if [ -e ${pathResults}/${motif} ]; then
    echo "ENCODE bg done"
else
   fimo -o ${pathResults}/${motif} ${pathJASPAR}/${motif}.meme ${pathResults}/ENCODE.Laverre2022.random50K.fa
fi

##########################################################
