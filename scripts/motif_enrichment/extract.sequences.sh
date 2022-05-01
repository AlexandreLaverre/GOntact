#!/bin/bash

export sp=$1
export genome=$2
export acc=$3
export dir=$4

##########################################################

pathResults=../../results/motif_enrichment/${sp}
pathGenome=../../data/genome_sequences/${sp}

##########################################################

perl extract.sequences.pl --pathGenomeSequence=${pathGenome}/${genome}.fa --pathCoordinates=${pathResults}/${acc}/${dir}/associated_enhancers.bed --coordConvention=0_open_end --pathOutput=${pathResults}/${acc}/${dir}/associated_enhancers.fa

##########################################################
## check if bg done

if [ -e ${pathResults}/ENCODE.Laverre2022.random50K.fa ]; then
    echo "ENCODE bg done"
else
    perl extract.sequences.pl --pathGenomeSequence=${pathGenome}/${genome}.fa --pathCoordinates=${pathResults}/ENCODE.Laverre2022.random50K.bed --coordConvention=0_open_end --pathOutput=${pathResults}/ENCODE.Laverre2022.random50K.fa
fi

##########################################################
