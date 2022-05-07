#!/bin/bash

#################################################################################

export sp=$1
export sample=$2
export domain=$3

#################################################################################

 ./_build/install/default/bin/gontact-utils compare-methods --path-annot1 <(zcat results/${sp}/${sample}/${domain}/GREAT_upstream5kb_downstream1kb_extend1Mb/GOntact_element_GO_association_background.txt.gz) --path-annot2 <(zcat results/${sp}/${sample}/${domain}/contacts_mindist0kb_maxdist1Mb/GOntact_element_GO_association_background.txt.gz) --path-output results/${sp}/${sample}/${domain}/compare_GO_assoc_background_GREAT_contacts_1Mb.txt

#################################################################################

 ./_build/install/default/bin/gontact-utils compare-methods --path-annot1 <(zcat results/${sp}/${sample}/${domain}/GREAT_upstream5kb_downstream1kb_extend1Mb/GOntact_element_GO_association_background.txt.gz) --path-annot2 <(zcat results/${sp}/${sample}/${domain}/simulated_contacts_mindist0kb_maxdist1Mb/GOntact_element_GO_association_background.txt.gz) --path-output results/${sp}/${sample}/${domain}/compare_GO_assoc_background_GREAT_simulated_contacts_1Mb.txt


#################################################################################
