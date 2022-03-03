#!/bin/bash

######################################################################

export pathHiC="../../data/PCHi-C"

######################################################################

cut -f 1-3 ${pathHiC}/human/Digest_hg38_HindIII_None.txt.baitmap > ${pathHiC}/human/hg38.baitmap

cut -f 1-3 ${pathHiC}/mouse/Digest_mm10_HindIII_None.txt.baitmap > ${pathHiC}/mouse/mm10.baitmap

######################################################################
