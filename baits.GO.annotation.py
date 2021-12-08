#!/usr/bin/env python3
# coding=utf-8

import argparse
from collections import defaultdict

parser = argparse.ArgumentParser()
parser.add_argument("species", help="species")
parser.add_argument("--UniqueGO", action="store_true",
                    help="Get unique list of GO term associated to genes for each bait (default=complete)")
args = parser.parse_args()

######################################################################

GenomeAssembly = "hg38" if args.species == "human" else "mm10"
Prefix = "unique" if args.UniqueGO else "complete"

path = "/home/laverre/GOntact"
Baits = path + "/data/PCHi-C/" + args.species + "/bait_coords_" + GenomeAssembly + ".txt"
EnsemblAnnotation = path + "/data/ensembl_annotations/" + args.species + "/GeneNames_Ensembl94.txt"
AnnotatedGenes = path + "/data/GeneOntology/tmp.genes.annotation.txt"
OutputFile = path + "/results/" + args.species + "/" + Prefix + ".GO.annotated.baits.txt"

######################################################################
## Correspondence between EnsemblID and GeneName
GeneID_dic = {}
with open(EnsemblAnnotation, 'r') as f:
    for gene in f.readlines()[1:]:
        gene = gene.strip("\n")
        gene = gene.split("\t")
        GeneID = gene[0]
        GeneName = gene[1]

        GeneID_dic[GeneID] = GeneName

######################################################################
## GO Annotations of each genes
GeneGO_dic = defaultdict(list)
with open(AnnotatedGenes, 'r') as f:
    for gene in f.readlines()[1:]:
        gene = gene.strip("\n")
        gene = gene.split("\t")
        GeneName = gene[0]
        GOTerm = gene[1]

        GeneGO_dic[GeneName].append(GOTerm)

######################################################################
## Annotations of each bait & writing Output
Output = open(OutputFile, "w")
Output.write("BaitID\tGeneID\tGeneName\tGOID\n")

with open(Baits, 'r') as f:
    for bait in f.readlines()[1:]:
        bait = bait.strip("\n")
        bait = bait.split("\t")
        BaitID = bait[0]
        BaitGenesID = bait[5].split(",")

        GenesNames = [GeneID_dic[GeneID] for GeneID in BaitGenesID]
        GenesGOTerm = [GO for GO in GeneGO_dic[GeneName] for GeneName in GenesNames]

        if args.unique:
            GenesGOTerm = list(set(GenesGOTerm))

        Output.write(BaitID + "\t" + ",".join(BaitGenesID) + "\t" + ",".join(GenesNames) + "\t" + ",".join(GenesGOTerm) + '\n')

Output.close()




