#!/usr/bin/env python3
# coding=utf-8

import argparse
from collections import defaultdict
import pandas
import itertools
import os
import re
import argparse

########################################################################################################################
### Arguments ###
parser = argparse.ArgumentParser()
parser.add_argument("species", help="species")
parser.add_argument("Enhancers", help="Input file of enhancers coordinates in BED format")

# Options
parser.add_argument("--minDistance", nargs="?", default=0, const=0, type=int,
                    help="minimum distance (in bp) between bait and contacted region (default = 0bp)")
parser.add_argument("--maxDistance", nargs="?", default=2000000, const=0, type=int,
                    help="maximum distance (in bp) between bait and contacted region (default = 2Mb)")
parser.add_argument("--KeepBaitedEnhancers", action="store_true",
                    help="Keep contacts where Enhancers are in baits (default = False)")
parser.add_argument("--KeepTransContact", action="store_true",
                    help="Keep inter-chromosome and intra >2Mb contacts (default = False)")
parser.add_argument("--KeepBaitBait", action="store_true", help="Keep bait-bait contacts (default = False)")

args = parser.parse_args()

########################################################################################################################
### Define path and files ###

GenomeAssembly = "hg38" if args.species == "human" else "mm10"
Prefix = os.path.basename(args.Enhancers).strip('.bed')

path = "/beegfs/data/necsulea/GOntact/"
Baits = path + "/data/PCHi-C/" + args.species + "/bait_coords_" + GenomeAssembly + ".txt"
Fragments = path + "/data/PCHi-C/" + args.species + "/frag_coords_" + GenomeAssembly + ".txt"
Contacts = path + "/data/PCHi-C/" + args.species + "/all_interactions_simplified.txt"

InputEnhancers = path + "data/enhancers/" + args.species + "/" + args.Enhancers
PathOutput = path + "/results/" + args.species + "/" + Prefix
if not os.path.exists(PathOutput):
    os.makedirs(PathOutput)

Baited = ".BaitedEnh" if args.KeepBaitedEnhancers else ""
Trans = ".Trans" if args.KeepTransContact else ""
Bait2Bait = ".bait2bait" if args.KeepBaitBait else ""

ForegroundOutput = PathOutput + "/foreground.contacts" + Baited + Trans + Bait2Bait + ".txt"
BackgroundOutput = PathOutput + "/background.contacts" + Baited + Trans + Bait2Bait + ".txt"

########################################################################################################################
### Overlap between input enhancers and restriction fragments ###

# Create dictionary of coordinates
def Coord_Dict(Coordinates, type):
    dic = defaultdict(list)
    N = 0
    with open(Coordinates, 'r') as f:
        # Check if header is present
        start = 1 if bool(re.search('start', f.readline())) is True else 0
        f.seek(0)

        for i in f.readlines()[start:]:
            i = i.strip("\n")
            i = i.split("\t")

            # coord = (start, end, ID)
            coord = (int(i[2]), int(i[3]), str(i[0])) if type == "background" else (int(i[1]), int(i[2]), str(i[3]))
            chr = str(i[1]) if type == "background" else str(i[0])

            dic[chr].append(coord)
            N += 1

    # Sorting coordinates for each chromosome
    for k in dic.keys():
        dic[k] = list(set(tuple(x) for x in dic[k]))
        dic[k].sort(key=lambda x: x[0])

    return dic, N


EnhancersDict, NbEnhancers = Coord_Dict(InputEnhancers, "foreground")
print("Found", NbEnhancers, "input enhancers.")
FragmentsDict, NbFragments = Coord_Dict(Fragments, "background")

#### Attribute input enhancers to restriction fragments ####
# Overlap Enhancers to Fragments
OverlapDict = defaultdict(list)
OverlapEnh = []
MissingEnhancers = []

for chr in EnhancersDict.keys():
    first_i = 0
    for Enhancer in EnhancersDict[chr]:
        start = Enhancer[0]
        end = Enhancer[1]
        EnhancerID = str(Enhancer[2])

        if chr in FragmentsDict.keys():
            # Initialization of first possible overlapping interest position
            i = first_i
            while i < len(FragmentsDict[chr]) and FragmentsDict[chr][i][1] < start:
                i += 1
            first_i = i

            # Adding all overlapping interest position to reference position
            while i < len(FragmentsDict[chr]) and FragmentsDict[chr][i][0] <= end:
                OverlapEnh.append(EnhancerID)
                OverlapDict[FragmentsDict[chr][i][2]].append(EnhancerID)
                i += 1

        if EnhancerID not in OverlapEnh:
            MissingEnhancers.append(EnhancerID)

if len(MissingEnhancers) != 0:
    print("Warning:", len(MissingEnhancers), "enhancer(s) do not overlap any restriction fragments")
    print(MissingEnhancers)

# Get all fragments containing input enhancers
SelectedFragments = OverlapDict.keys()
print("Found", len(SelectedFragments), "associated restriction fragments.")

########################################################################################################################
### Get contacts involving input enhancers ###

# Get all contacts
AllContacts = pandas.read_csv(Contacts, sep='\t')

# Drop contacts outside of specified distance range
AllContacts = AllContacts.drop(AllContacts[(AllContacts.Synteny == "cis") & (AllContacts.Distance > args.maxDistance)].index)
AllContacts = AllContacts.drop(AllContacts[(AllContacts.Synteny == "cis") & (AllContacts.Distance < args.minDistance)].index)

# Drop contacts according to specified filters
if not args.KeepBaitBait:
    AllContacts = AllContacts.drop(AllContacts[AllContacts.Type == "baited"].index)

if not args.KeepTransContact:
    AllContacts = AllContacts.drop(AllContacts[(AllContacts.Synteny == "trans")].index)

print("Found", len(AllContacts.index), "background contacts.")

# Get all interactions involving selected fragments
SelectedContacts = AllContacts.loc[AllContacts['FragmentID'].isin(SelectedFragments)]

########################################################################################################################
# Print some stats
FragInContacts = list(set(SelectedContacts['FragmentID'].tolist()))
EnhInContactedFrag = [OverlapDict[frag] for frag in FragInContacts]
EnhInContactedFrag = set(list(itertools.chain(*EnhInContactedFrag)))

print("Found", len(SelectedContacts.index), "foreground contacts with", len(FragInContacts), "selected fragments",
      "from", len(EnhInContactedFrag), "input enhancers:", round(len(EnhInContactedFrag)/NbEnhancers, 2)*100,
      "% enhancers are present in contacted fragments.")

if args.KeepBaitedEnhancers:
    BaitedEnhancersContacts = AllContacts.loc[AllContacts['BaitID'].isin(SelectedFragments)]
    SelectedContacts = pandas.concat([SelectedContacts, BaitedEnhancersContacts]).drop_duplicates().reset_index(drop=True)

    BaitInContacts = list(set(BaitedEnhancersContacts['BaitID'].tolist()))
    EnhInBait = [OverlapDict[frag] for frag in BaitInContacts]
    EnhInBait = set(list(itertools.chain(*EnhInBait)))

    print("Found", len(BaitedEnhancersContacts.index), "foreground contacts with", len(BaitInContacts), "baited fragments",
          "from", len(EnhInBait), "input enhancers:", round(len(EnhInBait)/NbEnhancers, 2)*100, "% enhancers are present in baits.")


print("Writing output...")
AllContacts.to_csv(BackgroundOutput, sep="\t", index=False)
SelectedContacts.to_csv(ForegroundOutput, sep="\t", index=False)

print("Done!")

##############################################################################
