#!/bin/bash

########################################################################

export path=/beegfs/data/necsulea/GOntact
export pathGO=${path}/data/GeneOntology
export pathScripts=${path}/scripts/process_GO_annotations

########################################################################

for sp in human mouse
do	  
perl ${pathScripts}/identify.redundant.categories.pl  --pathGOAnnotations=${pathGO}/${sp}.gene.annotation.txt --pathOutputSimplifiedAnnotations=${pathGO}/${sp}.simplified.gene.annotation.txt --pathOutputGOClusters=${pathGO}/${sp}.GO.clusters.txt
done 

########################################################################
