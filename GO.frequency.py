#!/usr/bin/env python3
# coding=utf-8

import argparse
import pandas
import itertools

parser = argparse.ArgumentParser()
parser.add_argument("species", help="species")
parser.add_argument("--UniqueGO", action="store_true",
                    help="Get unique list of GO term associated to genes for each bait (default=complete)")
args = parser.parse_args()

######################################################################

AnnotationType = "unique" if args.UniqueGO else "complete"

path = "/beegfs/data/necsulea/GOntact/"
AnnotatedBaits = path + "/results/" + args.species + "/" + AnnotationType + ".GO.annotated.baits.txt"
ForegroundContacts = path + "/results/" + args.species + "/foreground.contacts.txt"
BackgroundContacts = path + "/results/" + args.species + "/background.contacts.txt"

OutputFile = path + "/results/" + args.species + "/" + AnnotationType + ".GO.frequency.txt"

######################################################################
# Dictionary of GO annotated Baits
Baits2GO = {}
with open(AnnotatedBaits, 'r') as f:
    for bait in f.readlines()[1:]:
        bait = bait.strip("\n")
        bait = bait.split("\t")
        BaitID = bait[0]
        GOIDs = bait[3]

        Baits2GO[BaitID] = GOIDs.split(',')

print("Found", len(Baits2GO), "baits with GO terms")

######################################################################

def GOTerm_Frequency(InputContacts):
    # Get list of interested baits
    Contacts = pandas.read_csv(InputContacts, sep='\t', usecols=["BaitID"])
    SelectedBaits = list(itertools.chain(*Contacts.values.tolist()))

    N = len(SelectedBaits)
    print("Found", N , "contacts.")

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
Output.write("GOTerm\tForegroundFrequency\tBackgroundFrequency\n")

for GOTerm in BackgroundFrequency.keys():
    if GOTerm != '':
        FreqForeground = str(ForegroundFrequency[GOTerm]) if GOTerm in ForegroundFrequency.keys() else str(0)
        Output.write(str(GOTerm) + "\t" + FreqForeground + "\t" + str(BackgroundFrequency[GOTerm]) + "\n")

Output.close()
print("Done!")

######################################################################


