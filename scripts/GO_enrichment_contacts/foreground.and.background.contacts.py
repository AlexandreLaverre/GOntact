#!/usr/bin/env python3
# coding=utf-8

import argparse
from collections import defaultdict
import pandas
import itertools
import os
import re
import argparse
pandas.options.mode.chained_assignment = None

########################################################################################################################
### Arguments ###
parser = argparse.ArgumentParser()
parser.add_argument("species", help="species")
parser.add_argument("Enhancers", help="Input file of enhancers coordinates in BED format")

# Options
parser.add_argument("--BackgroundEnhancers", nargs="?", help="Input file of background enhancers coordinates in BED format")
parser.add_argument("--ExtendOverlap", nargs="?", default=0, const=0, type=int,
                    help="allow greater overlap (in bp) between enhancers and restriction fragments (default = 0bp)")
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
print("###### Selection of background and foreground contacts ######")
print("Running with following parameters:", args)

########################################################################################################################
### Define path and files ###

GenomeAssembly = "hg38" if args.species == "human" else "mm10"
Prefix = os.path.basename(args.Enhancers).strip('.bed')

path = "/beegfs/data/necsulea/GOntact/"  # "/home/laverre/Documents/GOntact/"
Baits = path + "/data/PCHi-C/" + args.species + "/bait_coords_" + GenomeAssembly + ".txt"
Fragments = path + "/data/PCHi-C/" + args.species + "/frag_coords_" + GenomeAssembly + ".txt"
Contacts = path + "/data/PCHi-C/" + args.species + "/all_interactions_simplified.txt"

InputEnhancers = path + "data/enhancers/" + args.species + "/" + args.Enhancers.strip('.bed') + ".bed"
if args.BackgroundEnhancers:
    BackgroundEnhancers = path + "data/enhancers/" + args.species + "/" + args.BackgroundEnhancers.strip('.bed') + ".bed"

minDistance = "mindist" + str(int(args.minDistance/1000)) + "Kb"
maxDistance = "_maxdist" + str(int(args.maxDistance/1000000)) + "Mb"
extendOverlap = "_extendOverlap" + str(int(args.ExtendOverlap/1000)) + "Kb"

PathOutput = path + "/results/GO_enrichment_contacts/" + args.species + "/" +\
             minDistance + maxDistance + extendOverlap + "/" + Prefix + "/ContactsSets"
if not os.path.exists(PathOutput):
    os.makedirs(PathOutput)

BackgroundType = args.BackgroundEnhancers.strip('.bed') if args.BackgroundEnhancers else "AllContacts"
Baited = ".BaitedEnh" if args.KeepBaitedEnhancers else ""
Trans = ".Trans" if args.KeepTransContact else ""
Bait2Bait = ".bait2bait" if args.KeepBaitBait else ""

ForegroundOutput = PathOutput + "/Foreground.Contacts" + Baited + Trans + Bait2Bait + ".txt"
BackgroundOutput = PathOutput + "/Background." + BackgroundType + Baited + Trans + Bait2Bait + ".txt"

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

            if type == "RestrictionFragments":
                chr, start, end = str(i[1]), int(i[2]), int(i[3])
            else:
                chr, start, end = str(i[0]), int(i[1]), int(i[2])

            ID = str(chr + ':' + str(start) + ':' + str(end))

            coord = (start, end, ID)
            dic[chr].append(coord)

            N += 1

    # Sorting coordinates for each chromosome
    for k in dic.keys():
        dic[k] = list(set(tuple(x) for x in dic[k]))
        dic[k].sort(key=lambda x: x[0])

    return dic, N


EnhancersDict, NbEnhancers = Coord_Dict(InputEnhancers, "ForegroundEnhancers")
print("Found", NbEnhancers, "input enhancers.")
FragmentsDict, NbFragments = Coord_Dict(Fragments, "RestrictionFragments")

if args.BackgroundEnhancers:
    BackgroundEnhancersDict, NbBackgroundEnhancers = Coord_Dict(BackgroundEnhancers, "BackgroundEnhancers")
    print("Found", NbBackgroundEnhancers, "background enhancers.")

#### Attribute input enhancers to restriction fragments ####
# Overlap Enhancers to Fragments

def OverlapEnhToFrag(EnhancersCoord, type):
    OverlapDict = defaultdict(list)
    OverlapEnh = []
    MissingEnhancers = []

    for chr in EnhancersCoord.keys():
        first_i = 0
        for Enhancer in EnhancersCoord[chr]:
            start = Enhancer[0]
            end = Enhancer[1]
            EnhancerID = str(Enhancer[2])
            notoverlap = True

            if chr in FragmentsDict.keys():
                # Initialization of first possible overlapping interest position
                i = first_i
                while i < len(FragmentsDict[chr]) and FragmentsDict[chr][i][1] < (start - args.ExtendOverlap):
                    i += 1
                first_i = i

                # Adding all overlapping interest position to reference position
                while i < len(FragmentsDict[chr]) and FragmentsDict[chr][i][0] <= (end + args.ExtendOverlap):
                    OverlapEnh.append(EnhancerID)
                    OverlapDict[FragmentsDict[chr][i][2]].append(EnhancerID)
                    notoverlap = False
                    i += 1

            if notoverlap:
                MissingEnhancers.append(EnhancerID)

    if len(MissingEnhancers) != 0:
        print("Warning:", len(MissingEnhancers), type, "enhancer(s) do not overlap any restriction fragments")
        print(MissingEnhancers)

    return OverlapDict


# Get all fragments containing input enhancers
InputOverlap = OverlapEnhToFrag(EnhancersDict, "input")
SelectedFragments = InputOverlap.keys()
print("Found", len(SelectedFragments), "restriction fragments associated with input enhancers.")

if args.BackgroundEnhancers:
    BackgroundOverlap = OverlapEnhToFrag(BackgroundEnhancersDict, "background")
    SelectedBackgroundFragments = BackgroundOverlap.keys()
    print("Found", len(SelectedBackgroundFragments), "restriction fragments associated with background enhancers.")

########################################################################################################################
### Get contacts involving input enhancers ###

# Get all contacts
AllContacts = pandas.read_csv(Contacts, sep='\t')

print("Found", len(AllContacts.index), "total contacts.")

# Drop contacts outside of specified distance range
AllContacts = AllContacts.drop(AllContacts[(AllContacts.Synteny == "cis") & (AllContacts.Distance > args.maxDistance)].index)
AllContacts = AllContacts.drop(AllContacts[(AllContacts.Synteny == "cis") & (AllContacts.Distance < args.minDistance)].index)

# Drop contacts according to specified filters
if not args.KeepBaitBait:
    AllContacts = AllContacts.drop(AllContacts[AllContacts.Type == "baited"].index)

if not args.KeepTransContact:
    AllContacts = AllContacts.drop(AllContacts[(AllContacts.Synteny == "trans")].index)

print("Found", len(AllContacts.index), "filtered contacts.")

# Get all interactions involving selected fragments
ForegroundContacts = AllContacts.loc[AllContacts['FragmentID'].isin(SelectedFragments)]

# Create new column for overlapped enhancers
ForegroundContacts['EnhancerID'] = ForegroundContacts['FragmentID'].map(InputOverlap).agg(','.join)

if args.BackgroundEnhancers:
    AllContacts = AllContacts.loc[AllContacts['FragmentID'].isin(SelectedBackgroundFragments)]
    AllContacts['EnhancerID'] = AllContacts['FragmentID'].map(BackgroundOverlap).agg(','.join)

########################################################################################################################
# Print some stats
ForegroundContactsFrag = list(set(ForegroundContacts['FragmentID'].tolist()))
ForegroundContactsEnh = [InputOverlap[frag] for frag in ForegroundContactsFrag]
ForegroundEnhinContact = set(list(itertools.chain(*ForegroundContactsEnh)))

print("Found", len(ForegroundContacts.index), "foreground contacts with", len(ForegroundContactsFrag), "selected fragments",
      "from", len(ForegroundEnhinContact), "input enhancers:", round((len(ForegroundEnhinContact)/NbEnhancers)*100, 2),
      "% enhancers are present in contacted fragments.")

if args.BackgroundEnhancers:
    BackgroundContactsFrag = list(set(AllContacts['FragmentID'].tolist()))
    BackgroundContactsEnh = [BackgroundOverlap[frag] for frag in BackgroundContactsFrag]
    BackgroundEnhinContact = set(list(itertools.chain(*BackgroundContactsEnh)))

    print("Found", len(AllContacts.index), "background contacts with", len(BackgroundContactsEnh), "selected fragments",
          "from", len(BackgroundEnhinContact), "background enhancers:", round((len(BackgroundEnhinContact) / NbBackgroundEnhancers) * 100, 2),
          "% enhancers are present in contacted fragments.")

if args.KeepBaitedEnhancers:
    BaitedEnhancersContacts = AllContacts.loc[AllContacts['BaitID'].isin(SelectedFragments)]
    SelectedContacts = pandas.concat([ForegroundContacts, BaitedEnhancersContacts]).drop_duplicates().reset_index(drop=True)

    BaitInContacts = list(set(BaitedEnhancersContacts['BaitID'].tolist()))
    EnhInBait = [InputOverlap[frag] for frag in BaitInContacts]
    EnhInBait = set(list(itertools.chain(*EnhInBait)))

    print("Found", len(BaitedEnhancersContacts.index), "foreground contacts with", len(BaitInContacts), "baited fragments",
          "from", len(EnhInBait), "input enhancers:", round((len(EnhInBait)/NbEnhancers)*100, 2), "% enhancers are present in baits.")


print("Writing output...")
AllContacts.to_csv(BackgroundOutput, sep="\t", index=False)
ForegroundContacts.to_csv(ForegroundOutput, sep="\t", index=False)

print("Contacts selection: Done!")

##############################################################################
