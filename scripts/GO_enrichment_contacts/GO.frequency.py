#!/usr/bin/env python3
# coding=utf-8

import argparse
import pandas
import itertools
import os

parser = argparse.ArgumentParser()
parser.add_argument("species", help="species")
parser.add_argument("Enhancers", help="species")
parser.add_argument("GONameSpace", help="GONameSpace : biological_process, molecular_function, cellular_component")

# Options
parser.add_argument("--ExtendOverlap", nargs="?", default=0, const=0, type=int,
                    help="allow greater overlap (in bp) between enhancers and restriction fragments (default = 0bp)")
parser.add_argument("--minDistance", nargs="?", default=0, const=0, type=int,
                    help="minimum distance (in bp) between bait and contacted region (default = 0bp)")
parser.add_argument("--maxDistance", nargs="?", default=2000000, const=0, type=int,
                    help="maximum distance (in bp) between bait and contacted region (default = 2Mb)")

parser.add_argument("--EnhancerContact", action="store_true", help="Get contacts at bait-enhancer scale (default = False)")
parser.add_argument("--EnhancerCountOnce", action="store_true", help="Enhancer can contribute only once for a given GOTerm (default = False)")
parser.add_argument("--BackgroundEnhancers", nargs="?", help="Input file of background enhancers coordinates in BED format")

parser.add_argument("--UniqueGO", action="store_true",
                    help="Get unique list of GO term associated to genes for each bait (default=complete)")
parser.add_argument("--KeepBaitedEnhancers", action="store_true", help="Keep contacts where Enhancers are in baits (default = False)")
parser.add_argument("--KeepTransContact", action="store_true", help="Keep inter-chromosome contacts (default = False)")
parser.add_argument("--KeepBaitBait", action="store_true", help="Keep inter-chromosome contacts (default = False)")

args = parser.parse_args()
print("### Estimate GOTerm frequency ###")
print("Running with following options:", args)

######################################################################
Prefix = os.path.basename(args.Enhancers).strip('.bed')

minDistance = "mindist" + str(int(args.minDistance/1000)) + "Kb"
maxDistance = "_maxdist" + str(int(args.maxDistance/1000000)) + "Mb"
extendOverlap = "_extendOverlap" + str(int(args.ExtendOverlap/1000)) + "Kb"

Prefix = minDistance + maxDistance + extendOverlap + "/" + Prefix

path = "/beegfs/data/necsulea/GOntact/"
pathResults = path + "/results/GO_enrichment_contacts/" + args.species
pathOutput = pathResults + "/" + Prefix + "/GOFrequency/"
if not os.path.exists(pathOutput):
    os.makedirs(pathOutput)

GODescription = path + "/data/GeneOntology/GODescription.txt"
AnnotationType = "complete." if not args.UniqueGO else ""
AnnotatedBaits = pathResults + "/" + AnnotationType + "GO.annotated.baits." + args.GONameSpace + ".txt"

# Define output files names according to options
EnhancerContact = ".ContactScale" if not args.EnhancerContact else ""
BackgroundType = args.BackgroundEnhancers.strip('.bed') if args.BackgroundEnhancers else "AllContacts"
EnhancerContribution = ".EnhAllCount" if not args.EnhancerCountOnce else ""
Baited = ".BaitedEnh" if args.KeepBaitedEnhancers else ""
Trans = ".Trans" if args.KeepTransContact else ""
Bait2Bait = ".bait2bait" if args.KeepBaitBait else ""

ForegroundContacts = pathResults + "/" + Prefix + "/ContactsSets/Foreground.Contacts" + Baited + Trans + Bait2Bait + ".txt"
BackgroundContacts = pathResults + "/" + Prefix + "/ContactsSets/Background." + BackgroundType + Baited + Trans + Bait2Bait + ".txt"

OutputFile = pathOutput + AnnotationType + args.GONameSpace + ".Background." + BackgroundType +\
             EnhancerContact + EnhancerContribution + Baited + Trans + Bait2Bait + ".txt"

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

        if GOIDs != "":
            Baits2GO[BaitID] = GOIDs.split(',')

print("Found", len(Baits2GO), "baits with GO terms")

######################################################################
# Calculate the frequency of each contacted GOTerm

def GOTerm_Frequency(InputContacts, type):
    scale = "enhancer" if (type == "foreground" and args.EnhancerContact) or (type == "background" and args.BackgroundEnhancers) else "fragment"

    # Get list of interested baits
    if type == "foreground" or args.BackgroundEnhancers:
        Contacts = pandas.read_csv(InputContacts, sep='\t', usecols=["BaitID", "EnhancerID"])
    else:
        Contacts = pandas.read_csv(InputContacts, sep='\t', usecols=["BaitID", "FragmentID"])

    print("Found", len(Contacts), " bait-fragment contacts.")

    # Count number of GO Term in InputContacts
    GOFrequency = {}
    NbContact = 0

    for contact in range(len(Contacts)):
        baitID = Contacts.iloc[contact, 0]
        ContactedID = Contacts.iloc[contact, 1]

        if baitID in Baits2GO.keys():  # Baits may don't have any associated GO Term
            NbContact += 1
            for GOID in Baits2GO[baitID]:
                GOFrequency.setdefault(GOID, [])

                # Bait-enhancer scale
                if scale == "enhancer":
                    for enhancer in ContactedID.split(','):
                        GOFrequency[GOID].append(enhancer)

                # Bait-fragment scale
                else:
                    GOFrequency[GOID].append(ContactedID)

    # Count total number of GOTerm-ContactedID contact
    if args.EnhancerCountOnce:
        for GOID in GOFrequency.keys():
            GOFrequency[GOID] = list(set(GOFrequency[GOID]))

    NbTotalContactedID = len(set(list(itertools.chain(*GOFrequency.values()))))

    print("Found", len(GOFrequency.keys()), "associated GOTerm.")
    print("Found", NbContact, " bait-other contacts with GOTerm.")
    print("Found", NbTotalContactedID, "contacted", scale, ".")

    return GOFrequency, NbTotalContactedID


print("## Running Foreground:")
ForegroundFrequency, NbForeground = GOTerm_Frequency(ForegroundContacts, "foreground")
print("## Running Background:")
BackgroundFrequency, NbBackground = GOTerm_Frequency(BackgroundContacts, "background")

######################################################################
print("Writing output...")
Output = open(OutputFile, "w")
Output.write("# Total number of Foreground contacts = " + str(NbForeground) + "\n")
Output.write("# Total number of Background contacts = " + str(NbBackground) + "\n")
Output.write("GOName\tGOTerm\tForegroundFrequency\tBackgroundFrequency\tForegroundEnhancers\n")

for GOTerm in BackgroundFrequency.keys():
    if GOTerm != '':
        GOName = GOInfos[GOTerm][0] if GOTerm in GOInfos else "Cluster"
        # Background Frequency
        FreqBackground = str(len(BackgroundFrequency[GOTerm]))

        # Foreground Frequency
        FreqForeground = str(0)
        EnhancersForeground = "NA"
        if GOTerm in ForegroundFrequency.keys():
            FreqForeground = str(len(ForegroundFrequency[GOTerm]))
            EnhancersForeground = ','.join(ForegroundFrequency[GOTerm])

        Output.write(GOName + "\t" + str(GOTerm) + "\t" + FreqForeground + "\t" + FreqBackground + "\t" + EnhancersForeground + "\n")

Output.close()
print("Done!")

######################################################################
