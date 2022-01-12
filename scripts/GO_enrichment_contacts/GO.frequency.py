#!/usr/bin/env python3
# coding=utf-8

import argparse
import pandas
import itertools
import os

parser = argparse.ArgumentParser()
parser.add_argument("species", help="species")
parser.add_argument("Enhancers", help="species")

# Options
parser.add_argument("--EnhancerContact", action="store_true", help="Get contacts at bait-enhancer scale (default = False)")
parser.add_argument("--BackgroundEnhancers", action="store_true", help="Consider the background set as annotated with enhancers (default = False)")
parser.add_argument("--UniqueGO", action="store_true",
                    help="Get unique list of GO term associated to genes for each bait (default=complete)")
parser.add_argument("--EnhancerCountOnce", action="store_true", help="Enhancer can contribute only once for a given GOTerm (default = False)")
parser.add_argument("--WithoutPropagation", action="store_true", help="Get gene annotation without GO propagation")
parser.add_argument("--KeepBaitedEnhancers", action="store_true", help="Keep contacts where Enhancers are in baits (default = False)")
parser.add_argument("--KeepTransContact", action="store_true", help="Keep inter-chromosome contacts (default = False)")
parser.add_argument("--KeepBaitBait", action="store_true", help="Keep inter-chromosome contacts (default = False)")

args = parser.parse_args()

######################################################################
Prefix = os.path.basename(args.Enhancers).strip('.bed')
AnnotationType = "unique" if args.UniqueGO else "complete"

path = "/home/laverre/Documents/GOntact/" #"/beegfs/data/necsulea/GOntact/"
pathResults = path + "/results/GO_enrichment_contacts/" + args.species

GODescription = path + "/data/GeneOntology/GODescription.txt"

if args.WithoutPropagation:
    AnnotatedBaits = pathResults + "/before.propagation/" + AnnotationType + ".GO.annotated.baits.txt"
else:
    AnnotatedBaits = pathResults + "/" + AnnotationType + ".GO.annotated.baits.txt"

# Define output files names according to options
EnhancerContact = ".EnhancerScale" if args.EnhancerContact else ""
BackgroundType = "Background.EnhancersContacts" if args.BackgroundEnhancers else "Background.AllContacts"
EnhancerContribution = ".EnhancerCountOnce" if args.EnhancerCountOnce else ""
Baited = ".BaitedEnh" if args.KeepBaitedEnhancers else ""
Trans = ".Trans" if args.KeepTransContact else ""
Bait2Bait = ".bait2bait" if args.KeepBaitBait else ""
WithoutPropagation = ".WithoutPropagation" if args.WithoutPropagation else ""

ForegroundContacts = pathResults + "/" + Prefix + "/Foreground.Contacts" + Baited + Trans + Bait2Bait + ".txt"
BackgroundContacts = pathResults + "/" + Prefix + "/" + BackgroundType + Baited + Trans + Bait2Bait + ".txt"

OutputFile = pathResults + "/" + Prefix + "/" + AnnotationType + ".GO.frequency." + BackgroundType + EnhancerContact \
             + EnhancerContribution + Baited + Trans + Bait2Bait + WithoutPropagation + ".txt"

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

def GOTerm_Frequency(InputContacts, type):
    # Get list of interested baits
    if type == "foreground" or args.BackgroundEnhancers:
        Contacts = pandas.read_csv(InputContacts, sep='\t', usecols=["BaitID", "EnhancerID"])
    else:
        Contacts = pandas.read_csv(InputContacts, sep='\t', usecols=["BaitID", "FragmentID"])

    N = len(Contacts)
    print("Found", N, " bait-fragment contacts.")

    # List of baits
    SelectedBaits = Contacts["BaitID"].tolist()

    # Count number of GO Term in InputContacts
    GOFrequency = {}
    for BaitID in SelectedBaits:
        if BaitID in Baits2GO.keys():           # Baits may don't have any associated GO Term
            for GOID in Baits2GO[BaitID]:
                GOFrequency.setdefault(GOID, 0)
                GOFrequency[GOID] += 1

    print("Found", len(GOFrequency.keys()), "associated GOTerm.")

    # Contacts between baits and enhancers
    if (type == "foreground" and args.EnhancerContact) or (type == "background" and args.BackgroundEnhancers):
        ContactByEnhancer = {}
        N = 0
        for x in range(len(Contacts)):
            baitID = Contacts.iloc[x, 0]
            enhancersID = Contacts.iloc[x, 1]

            if enhancersID != 'NA':
                ContactByEnhancer.setdefault(baitID, [])
                for enhancerID in enhancersID.split(','):
                    if enhancerID not in ContactByEnhancer[baitID]:
                        ContactByEnhancer[baitID].append(enhancerID)
                        N += 1

        print("Found", N, " bait-enhancer contacts.")

        GOFrequency = {}
        for BaitID in ContactByEnhancer.keys():
            if BaitID in Baits2GO.keys():  # Baits may don't have any associated GO Term
                for GOID in Baits2GO[BaitID]:
                    GOFrequency.setdefault(GOID, [])
                    for enhancerID in ContactByEnhancer[BaitID]:
                        # Enhancer contribution can be binary or not
                        if args.EnhancerCountOnce and enhancerID not in GOFrequency[GOID]:
                            GOFrequency[GOID].append(enhancerID)
                        else:
                            GOFrequency[GOID].append(enhancerID)

    return GOFrequency, N


print("## Running Foreground:")
ForegroundFrequency, NbForeground = GOTerm_Frequency(ForegroundContacts, "foreground")
print("## Running Background:")
BackgroundFrequency, NbBackground = GOTerm_Frequency(BackgroundContacts, "background")

######################################################################
print("Writing output...")
Output = open(OutputFile, "w")
Output.write("# Total number of Foreground contacts = " + str(NbForeground) + "\n")
Output.write("# Total number of Background contacts = " + str(NbBackground) + "\n")
Output.write("GOName\tGONamespace\tGOTerm\tForegroundFrequency\tBackgroundFrequency\n")

for GOTerm in BackgroundFrequency.keys():
    if GOTerm != '':
        FreqForeground = str(0)
        if GOTerm in ForegroundFrequency.keys():
            if args.EnhancerContact:
                ForegroundFrequency[GOTerm] = len(ForegroundFrequency[GOTerm])
            FreqForeground = str(ForegroundFrequency[GOTerm])

        if args.BackgroundEnhancers:
            BackgroundFrequency[GOTerm] = len(BackgroundFrequency[GOTerm])
        FreqBackground = str(BackgroundFrequency[GOTerm])

        Output.write(GOInfos[GOTerm][0] + "\t" + GOInfos[GOTerm][1] + "\t" + str(GOTerm) + "\t" +
                     FreqForeground + "\t" + FreqBackground + "\n")

Output.close()
print("Done!")

######################################################################
