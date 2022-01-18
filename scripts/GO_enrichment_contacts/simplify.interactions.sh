#!/bin/bash

export species="mouse"	# i.e: human or mouse

export path="/beegfs/data/necsulea/GOntact/data/PCHi-C/${species}/"

######################################################

# Format ID
cut -f 1-3 ${path}/all_interactions.txt > ${path}/BaitID
cut -f 4-6 ${path}/all_interactions.txt > ${path}/FragmentID

sed -i "s/\t/:/g" ${path}/BaitID
sed -i "s/\t/:/g" ${path}/FragmentID

# Get Contact Type and distance
cut -f 7-8 ${path}/all_interactions.txt > ${path}/Stats

# Get the Synteny
awk '{if ($1 == $4) {print "cis"} else {print "trans"} }' ${path}/all_interactions.txt > ${path}/Synteny

# Concatenate infos and change header
paste -d "\t" ${path}/BaitID ${path}/FragmentID ${path}/Stats ${path}/Synteny > ${path}/all_interactions_simplified.txt
sed -i "1s/.*/BaitID\tFragmentID\tType\tDistance\tSynteny/" ${path}/all_interactions_simplified.txt

rm ${path}/BaitID ${path}/FragmentID ${path}/Stats ${path}/Synteny

######################################################
