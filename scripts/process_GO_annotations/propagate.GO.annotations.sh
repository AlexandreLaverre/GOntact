#!/bin/bash

########################################################################

export path=/beegfs/data/necsulea/GOntact
export pathGO=${path}/data/GeneOntology
export pathScripts=${path}/scripts/process_GO_annotations

########################################################################
## human

perl ${pathScripts}/propagate.GO.annotations.pl --pathGOTree=${pathGO}/go-basic.obo --pathInputGeneAnnotation=${pathGO}/goa_human.gaf  --pathOutputGeneAnnotation=${pathGO}/human.gene.annotation.txt --pathOutputGOCategories=${pathGO}/GOCategories.txt

########################################################################
## mouse

perl ${pathScripts}/propagate.GO.annotations.pl --pathGOTree=${pathGO}/go-basic.obo --pathInputGeneAnnotation=${pathGO}/mgi.gaf  --pathOutputGeneAnnotation=${pathGO}/mouse.gene.annotation.txt --pathOutputGOCategories=${pathGO}/GOCategories.txt

########################################################################
