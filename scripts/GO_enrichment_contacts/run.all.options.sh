#!/bin/bash
#Contacts options
export ExtendOverlap_options=(0 1000 2500)
export maxDistance_options=(1000000 2000000)
export minDistance_options=(0 25000)
export UniqueGO_options=(--UniqueGO "")
export BackgroundEnhancers_options=("" FANTOM5.selection.bed FANTOM5.Laverre2022.bed)

#GOFrequency options
export GONameSpace_options=(biological_process molecular_function cellular_component)
export EnhancerCountOnce_options=(--EnhancerCountOnce "")
export EnhancerContact_options=(--EnhancerContact "")
export Background_options=(--BackgroundEnhancers "")

#Count
numRound=0

for ExtendOverlap in "${ExtendOverlap_options[@]}"; do
for maxDistance in "${maxDistance_options[@]}"; do
for minDistance in "${minDistance_options[@]}"; do
for unique in "${UniqueGO_options[@]}"; do
for BackgroundEnhancer in "${BackgroundEnhancers_options[@]}"; do
./foreground.and.background.contacts.py human FANTOM5.first1000.kidney.enhancers.hg38.bed --ExtendOverlap ${ExtendOverlap} --maxDistance ${maxDistance} --minDistance ${minDistance} ${UniqueGO} --BackgroundEnhancers ${BackgroundEnhancer}

#GOFrequency
for GOName in "${GONameSpace_options[@]}"; do
for EnhancerCountOnce in "${EnhancerCountOnce_options[@]}"; do
for EnhancerContact in "${EnhancerContact_options[@]}"; do
./GO.frequency.py human FANTOM5.first1000.kidney.enhancers.hg38.bed ${GOName} --ExtendOverlap ${ExtendOverlap} --maxDistance ${maxDistance} --minDistance ${minDistance} ${UniqueGO} ${EnhancerCountOnce} ${EnhancerContact}

if [[ ${BackgroundEnhancer} != "" ]] & [[ ${EnhancerContact} == "--EnhancerContact" ]]; then
./GO.frequency.py human FANTOM5.first1000.kidney.enhancers.hg38.bed ${GOName} --ExtendOverlap ${ExtendOverlap} --maxDistance ${maxDistance} --minDistance ${minDistance} ${UniqueGO} ${EnhancerCountOnce} ${EnhancerContact} --BackgroundEnhancers ${BackgroundEnhancer}
fi
numRound=$((numRound+1))
echo "Number of Round: $numRound"
echo "##########################"

done
done
done
done
done
done
done
done


