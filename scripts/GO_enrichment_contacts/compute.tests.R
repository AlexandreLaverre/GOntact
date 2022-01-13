#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

library(data.table)

path="/beegfs/data/necsulea/GOntact/results/GO_enrichment_contacts/human/"
dist = args[1]
file = args[2]

GOFile = paste0(path, dist, "/FANTOM5.first1000.kidney.enhancers.hg38/GOFrequency/", file)

## Get the total Foreground & Background count in the header
HeadFile <- file(GOFile,"r")
ForegroundCount = as.numeric(sub(".*= ", "",  readLines(HeadFile,n=1)))
BackgroundCount = as.numeric(sub(".*= ", "",  readLines(HeadFile,n=1)))
close(HeadFile)

GOEnrichment = fread(GOFile, h=T, quote="", sep="\t", select = c(1:4))

# Hypergeometric test
GOEnrichment$HyperPval = apply(GOEnrichment, 1, function(x) phyper(as.numeric(x["ForegroundFrequency"])-1, as.numeric(x["BackgroundFrequency"]), 
                                                              BackgroundCount-as.numeric(x["BackgroundFrequency"]), ForegroundCount, lower.tail=FALSE))

#GOEnrichment$BinomPval = apply(GOEnrichment, 1, function(x) phyper(as.numeric(x["ForegroundFrequency"])-1, as.numeric(x["BackgroundFrequency"]), 
#                                                                   BackgroundCount-as.numeric(x["BackgroundFrequency"]), ForegroundCount, lower.tail=FALSE))

# FDR                   
GOEnrichment$FDR = p.adjust(GOEnrichment$HyperPval, method = "fdr")

# Output

GOEnrichment <- GOEnrichment[order(GOEnrichment$FDR),]
write.table(GOEnrichment, file = paste0(path, dist, "/FANTOM5.first1000.kidney.enhancers.hg38/GOEnrichment/", file),
            quote=FALSE, sep = '\t', dec = '.', row.names = F, col.names = T)


# #######################################################################################################
# # Distance distribution between fragments
# 
# par(mfrow=c(2,2))
# for (sample in samples){
#   enh=paste("FANTOM5.first1000", sample, "enhancers.hg38", sep=".")
#   
#   ForegroundContact = read.table(paste0(path, sp, "/", enh, "/foreground.contacts.txt"), h=T, sep="\t")
#   hist(ForegroundContact$Distance, breaks=50, main=sample, xlab="Distance")
#   med = median(ForegroundContact$Distance)
#   abline(v=med, col="red")
#   mtext(paste0("med.dist=", med, "\n", "N=", nrow(ForegroundContact)), line = -1)
#   
# }