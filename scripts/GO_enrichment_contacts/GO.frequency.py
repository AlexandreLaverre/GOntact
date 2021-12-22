#!/usr/bin/env python3
# coding=utf-8

import argparse
import pandas
import itertools
import os

parser = argparse.ArgumentParser()
parser.add_argument("species", help="species")
parser.add_argument("Enhancers", help="species")
parser.add_argument("--WithoutPropagation", action="store_true", help="Get gene annotation without GO propagation")
parser.add_argument("--KeepBaitedEnhancers", action="store_true", help="Keep contacts where Enhancers are in baits (default = False)")
parser.add_argument("--KeepTransContact", action="store_true", help="Keep inter-chromosome contacts (default = False)")
parser.add_argument("--KeepBaitBait", action="store_true", help="Keep inter-chromosome contacts (default = False)")
parser.add_argument("--UniqueGO", action="store_true",
                    help="Get unique list of GO term associated to genes for each bait (default=complete)")
parser.add_argument("--CountBaitOnce", action="store_true",
                    help="Count only one time each bait on foreground and background set")
args = parser.parse_args()

######################################################################
Prefix = os.path.basename(args.Enhancers).strip('.bed')
AnnotationType = "unique" if args.UniqueGO else "complete"

path = "/beegfs/data/necsulea/GOntact/"
pathResults = path + "/results/" + args.species

GODescription = path + "/data/GeneOntology/GODescription.txt"

if args.WithoutPropagation:
    AnnotatedBaits = pathResults + "/before.propagation/" + AnnotationType + ".GO.annotated.baits.txt"
else:
    AnnotatedBaits = pathResults + "/" + AnnotationType + ".GO.annotated.baits.txt"

Baited = ".BaitedEnh" if args.KeepBaitedEnhancers else ""
Trans = ".Trans" if args.KeepTransContact else ""
Bait2Bait = ".bait2bait" if args.KeepBaitBait else ""
CountOnce = ".CountOnce" if args.CountBaitOnce else ""
WithoutPropagation = ".WithoutPropagation" if args.WithoutPropagation else ""

ForegroundContacts = pathResults + "/" + Prefix + "/foreground.contacts" + Baited + Trans + Bait2Bait + ".txt"
BackgroundContacts = pathResults + "/" + Prefix + "/background.contacts" + Baited + Trans + Bait2Bait + ".txt"

OutputFile = pathResults + "/" + Prefix + "/" + AnnotationType + ".GO.frequency" + Baited + Trans + \
             Bait2Bait + CountOnce + WithoutPropagation + ".txt"

######################################################################
# Dictionary of GO Description
GOInfos = {}
with open(GODescription, 'r') as f:
    for GO in f.readlines()[1:]:
        GO = GO.strip("\n")
        GO = GO.split("\t")
        GOID, GOName, GONamespace = GO[0], GO[1], GO[2]

        GOInfos[GOID] = (GOName, GONamespace)

######################################################################
# Dictionary of GO annotated Baits
Baits2GO = {}
with open(AnnotatedBaits, 'r') as f:
    for bait in f.readlines()[1:]:
        bait = bait.strip("\n")
        bait = bait.split("\t")
        BaitID, GOIDs = bait[0], bait[3]

        Baits2GO[BaitID] = GOIDs.split(',')

print("Found", len(Baits2GO), "baits with GO terms")

######################################################################

def GOTerm_Frequency(InputContacts):
    # Get list of interested baits
    Contacts = pandas.read_csv(InputContacts, sep='\t', usecols=["BaitID"])
    SelectedBaits = list(itertools.chain(*Contacts.values.tolist()))

    if args.CountBaitOnce:
        SelectedBaits = list(set(SelectedBaits))

    N = len(SelectedBaits)
    print("Found", N, "contacts.")

    # Count number of GO Term in InputContacts
    GOFrequency = {}
    for BaitID in SelectedBaits:
        if BaitID in Baits2GO.keys():           # Baits may don't have any associated GO Term
            for GOID in Baits2GO[BaitID]:
                if GOID in GOFrequency.keys():
                    GOFrequency[GOID] = GOFrequency[GOID] + 1
                else:
                    GOFrequency[GOID] = 1

    print("Found", len(GOFrequency.keys()), "associated GOTerm.")

    return GOFrequency, N


print("## Running Foreground:")
ForegroundFrequency, NbForeground = GOTerm_Frequency(ForegroundContacts)
print("## Running Background:")
BackgroundFrequency, NbBackground = GOTerm_Frequency(BackgroundContacts)

######################################################################
print("Writing output...")
Output = open(OutputFile, "w")
Output.write("# Total number of Foreground contacts = " + str(NbForeground) + "\n")
Output.write("# Total number of Background contacts = " + str(NbBackground) + "\n")
Output.write("GOName\tGONamespace\tGOTerm\tForegroundFrequency\tBackgroundFrequency\n")

for GOTerm in BackgroundFrequency.keys():
    if GOTerm != '':
        FreqForeground = str(ForegroundFrequency[GOTerm]) if GOTerm in ForegroundFrequency.keys() else str(0)
        Output.write(GOInfos[GOTerm][0] + "\t" + GOInfos[GOTerm][1] + "\t" + str(GOTerm) + "\t" +
                     FreqForeground + "\t" + str(BackgroundFrequency[GOTerm]) + "\n")

Output.close()
print("Done!")

######################################################################


