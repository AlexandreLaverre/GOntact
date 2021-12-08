#!/usr/bin/env python3
# coding=utf-8

import argparse
import pandas
import itertools

parser = argparse.ArgumentParser()
parser.add_argument("species", help="species")
parser.add_argument("ContactType", help="foreground or background")
parser.add_argument("--UniqueGO", action="store_true",
                    help="Get unique list of GO term associated to genes for each bait (default=complete)")
args = parser.parse_args()

######################################################################

AnnotationType = "unique" if args.UniqueGO else "complete"

path = "/beegfs/data/necsulea/GOntact/"
AnnotatedBaits = path + "/results/" + args.species + "/" + AnnotationType + ".GO.annotated.baits.txt"
InputContacts = path + "/results/" + args.species + "/" + args.ContactType + ".contacts.txt"

OutputFile = path + "/results/" + args.species + "/" + args.ContactType + "." + AnnotationType + ".GO.frequency.txt"

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
# Get list of interested baits
Contacts = pandas.read_csv(InputContacts, sep='\t', usecols=["BaitID"])
SelectedBaits = list(itertools.chain(*Contacts.values.tolist()))

print("Found", len(SelectedBaits), "contacts.")

# Count number of GO Term in InputContacts
GOFrequency = {}
for BaitID in SelectedBaits:
    if BaitID in Baits2GO.keys():           # Daits may don't have any associated GO Term
        for GOID in Baits2GO[BaitID]:
            if GOID in GOFrequency.keys():
                GOFrequency[GOID] = GOFrequency[GOID] + 1
            else:
                GOFrequency[GOID] = 1

print("Found", len(GOFrequency.keys()), "associated GOTerm.")

######################################################################
print("Writing output...")
Output = open(OutputFile, "w")
Output.write("# Total number of contacts = " + str(len(SelectedBaits)) + "\n")
Output.write("GOTerm\tFrequency\n")

for GOTerm, Frequency in GOFrequency.items():
    if GOTerm != '':
        Output.write(str(GOTerm) + "\t" + str(Frequency) + "\n")

Output.close()
print("Done!")

######################################################################


