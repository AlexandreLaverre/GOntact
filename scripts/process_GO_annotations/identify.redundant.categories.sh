#!/bin/bash

########################################################################

export path=/beegfs/data/necsulea/GOntact
export pathGO=${path}/data/GeneOntology
export pathScripts=${path}/scripts/process_GO_annotations

########################################################################

for sp in human mouse
do
    for space in biological_process molecular_function cellular_component
    do
	perl ${pathScripts}/identify.redundant.categories.pl  --pathGOCategories=${pathGO}/GOCategories.txt --GOSpace=${space} --pathGOAnnotations=${pathGO}/${sp}.gene.annotation.txt --pathOutputSimplifiedAnnotations=${pathGO}/${sp}.simplified.gene.annotation.${space}.txt --pathOutputGOClusters=${pathGO}/${sp}.GO.clusters.${space}.txt
    done
done

########################################################################
