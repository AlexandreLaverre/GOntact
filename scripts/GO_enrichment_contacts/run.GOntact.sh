#!/bin/bash

export species=$1			# i.e: human or mouse
export InputEnhancers=$2		# i.e: FANTOM5.first1000.kidney.enhancers.hg38
export BackgroundEnhancers=$3		# i.e: FANTOM5.Laverre2022 or FANTOM5.selection
export minDist=$4			# i.e: 0 or 25000
export maxDist=$5			# i.e: 1000000 or 2000000
export ExtendOverlap=$6			# i.e: 0 or 5000

export GONameSpace="biological_process"

######################################################
minDistname=$(($minDist/1000))
maxDistname=$(($maxDist/1000000))
overlapname=$(($ExtendOverlap/1000))
foldername="mindist${minDistname}Kb_maxdist${maxDistname}Mb_extendOverlap${overlapname}Kb"

path="/beegfs/data/necsulea/GOntact/results/GO_enrichment_contacts/${species}"
pathResults="${path}/${foldername}/${InputEnhancers}"

pathScripts="/beegfs/data/alaverre/GOntact/scripts/GO_enrichment_contacts"
pythonPath="/beegfs/data/alaverre/Tools/envs/bin/python"

######################################################
### GO annotation of baits 
# Set to unique GOTerm for each bait (--UniqueGO option)

if [ -f "${path}/GO.annotated.baits.${GONameSpace}.txt" ]; then
    	echo "Baits already annotated."
else
	${pythonPath} ${pathScripts}/baits.GO.annotation.py ${species} ${GONameSpace} --UniqueGO
fi

######################################################
### Estimate Foreground and Background contacts 
# Set to avoid trans (--KeepTransContact) and bait-bait contacts (--KeepBaitBait), and don't consider baited enhancers (--KeepBaitedEnhancers))

if [ -f "${pathResults}/ContactsSets/Background.${BackgroundEnhancers}.txt" ]; then
    	echo "Foreground and Background contacts already determined."
else
	${pythonPath} ${pathScripts}/foreground.and.background.contacts.py ${species} ${InputEnhancers} --BackgroundEnhancers ${BackgroundEnhancers} --minDistance ${minDist} --maxDistance ${maxDist} --ExtendOverlap ${ExtendOverlap}
fi

######################################################
### GO Term frequency
# Set to Enhancer Scale (--EnhancerContact) and count enhancer only once for each GO Term (--EnhancerCountOnce)

if [ -f "${pathResults}/GOFrequency/${GONameSpace}.Background.${BackgroundEnhancers}.txt" ]; then
    	echo "GOFrequency already done."
else
	${pythonPath} ${pathScripts}/GO.frequency.py ${species} ${InputEnhancers} ${GONameSpace} --BackgroundEnhancers ${BackgroundEnhancers} --minDistance ${minDist} --maxDistance ${maxDist} --ExtendOverlap ${ExtendOverlap} --UniqueGO --EnhancerContact --EnhancerCountOnce 
fi

######################################################
### Compute Hypergeometric and Binomial tests
# Set to ordering output by Hypergeometric FDR

if [ -f "${pathResults}/GOEnrichment/${GONameSpace}.Background.${BackgroundEnhancers}.txt" ]; then
    	echo "Tests already done."
else
	mkdir -p "${pathResults}/GOEnrichment/"

	/bin/Rscript compute.tests.R  ${species} ${InputEnhancers} ${foldername} ${GONameSpace}.Background.${BackgroundEnhancers}.txt
fi

echo "All done!"

######################################################
