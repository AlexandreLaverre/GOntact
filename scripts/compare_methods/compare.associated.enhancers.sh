#!/bin/bash

export sp=$1
export dataset=$2
export method1=$3
export submethod1=$4
export method2=$5
export submethod2=$6
export background=$7
export space=$8

####################################################################################

export path=/beegfs/data/necsulea/GOntact
export pathGREAT=${path}/results/GO_enrichment_GREAT/${sp}
export pathGOntact=${path}/results/GO_enrichment_contacts/${sp}
export pathResults=${path}/results/method_comparison/${sp}/${dataset}
export pathScripts=${path}/scripts/compare_methods

export ensrelease=94

####################################################################################

if [ -e ${pathResults} ]; then
    echo "path output exists"
else
    mkdir -p ${pathResults}
fi

####################################################################################

if [ ${method1} = "GREAT" ]; then
    export path1=${pathGREAT}/${submethod1}/${dataset}/enhancer_GO_association_Ensembl${ensrelease}_${space}.txt 
fi

if [ ${method1} = "GOntact" ]; then
    export path1=${pathGOntact}/${submethod1}/${dataset}/GOFrequency/${space}.Background.${background}.txt
fi

####################################################################################

if [ ${method2} = "GREAT" ]; then
    export path2=${pathGREAT}/${submethod2}/${dataset}/enhancer_GO_association_Ensembl${ensrelease}_${space}.txt 
fi

if [ ${method2} = "GOntact" ]; then
    export path2=${pathGOntact}/${submethod2}/${dataset}/GOFrequency/${space}.Background.${background}.txt
fi

####################################################################################

perl ${pathScripts}/compare.associated.enhancers.pl --pathGOEnhancers1=${path1} --pathGOEnhancers2=${path2} --pathOutput=${pathResults}/${space}_${method1}.${submethod1}_vs_${method2}.${submethod2}.txt

####################################################################################
