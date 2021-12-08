#!/usr/bin/env python3
# coding=utf-8

import argparse
from collections import defaultdict
import pandas
import itertools

parser = argparse.ArgumentParser()
parser.add_argument("species", help="species")
parser.add_argument("InputEnhancers", help="Input file of enhancers coordinates in BED format")
parser.add_argument("--BaitedEnhancers", action="store_true",
                    help="Keep contacts where Enhancers are in baits (default = False)")
parser.add_argument("--TransContact", action="store_true",
                    help="Keep inter-chromosome contacts (default = False)")
args = parser.parse_args()

################################################

GenomeAssembly = "hg38" if args.species == "human" else "mm10"

path = "/home/laverre/GOntact"
Baits = path + "/data/PCHi-C/" + args.species + "/bait_coords_" + GenomeAssembly + ".txt"
Fragments = path + "/data/PCHi-C/" + args.species + "/frag_coords_" + GenomeAssembly + ".txt"
Contacts = path + "/data/PCHi-C/" + args.species + "/all_interactions.txt"
#OutputFile = path + "/results/" + args.species + "/foreground.contacts.txt"
OutputFile = "/home/laverre/Documents/GOntact/foreground.contacts.txt"

##############################################################################
######## Overlap between input enhancers and restriction fragments ########

# Create dictionary of coordinates
def Coord_Dict(Coordinates):
    dic = defaultdict(list)
    N = 0
    with open(Coordinates, 'r') as f:
        for i in f.readlines()[1:]:
            i = i.strip("\n")
            i = i.split("\t")
            chr = str(i[1])
            coord = (int(i[2]), int(i[3]), str(i[0]))  # coord = (start, end, ID)

            dic[chr].append(coord)
            N += 1

    # Sorting coordinates for each chromosome
    for k in dic.keys():
        dic[k] = list(set(tuple(x) for x in dic[k]))
        dic[k].sort(key=lambda x: x[0])

    return dic, N


FragmentsDict, NbFragments = Coord_Dict(Fragments)
EnhancersDict, NbEnhancers = Coord_Dict(args.InputEnhancers)
print("Found", NbEnhancers, "input enhancers.")

print("Attribute input enhancers to restriction fragments...")
# Overlap Enhancers to Fragments
OverlapDict = defaultdict(list)
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
                OverlapDict[EnhancerID].append(FragmentsDict[chr][i][2])
                i += 1

        if EnhancerID not in OverlapDict.keys():
            MissingEnhancers.append(EnhancerID)

if len(MissingEnhancers) != 0:
    print("Warning:", len(MissingEnhancers), "enhancers do not overlap any restriction fragments")

# Get all fragments containing input enhancers
SelectedFragments = list(set(list(itertools.chain(*OverlapDict.values()))))
print("Found", len(SelectedFragments), "associated restriction fragments.")

##############################################################################
print("Get contacts involving input enhancers...")

# Get all contacts
AllContacts = pandas.read_csv(Contacts, sep='\t', usecols=[0, 1, 2, 3, 4, 5, 7])

# Define IDs
AllContacts['BaitID'] = AllContacts[AllContacts.columns[0:3]].apply(lambda x: ':'.join(x.astype(str)), axis=1)
AllContacts['FragmentID'] = AllContacts[AllContacts.columns[3:6]].apply(lambda x: ':'.join(x.astype(str)), axis=1)

# Get all interactions involving selected fragments
if args.BaitedEnhancers:
    SelectedContacts = AllContacts.loc[
        AllContacts['FragmentID'].isin(SelectedFragments) | AllContacts['BaitID'].isin(SelectedFragments)]
else:
    SelectedContacts = AllContacts.loc[AllContacts['FragmentID'].isin(SelectedFragments)]

print("Found", len(SelectedContacts.index), "contacts.")

if not args.TransContact:
    SelectedContacts = SelectedContacts.loc[SelectedContacts['chr_bait'] == SelectedContacts['chr']
                                            & SelectedContacts['distance'] < 2000000]
    print("Found", len(SelectedContacts.index), "intra-chromosome & <2Mb contacts.")

print("Writing output...")
SelectedContacts.to_csv(OutputFile, sep="\t", index=False)

print("Done!")
