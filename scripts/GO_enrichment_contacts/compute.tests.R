##################################################################################
#!/usr/bin/env Rscript
library(data.table)

args = commandArgs(trailingOnly=TRUE)

dist = args[1]
file = args[2]

pathResult = "/beegfs/data/necsulea/GOntact/results/GO_enrichment_contacts/human/"
GOFile = paste0(pathResult, dist, "/FANTOM5.first1000.kidney.enhancers.hg38/GOFrequency/", file)

##################################################################################
## Get the total Foreground & Background count in the header
HeadFile <- file(GOFile,"r")
ForegroundCount = as.numeric(sub(".*= ", "",  readLines(HeadFile,n=1)))
BackgroundCount = as.numeric(sub(".*= ", "",  readLines(HeadFile,n=1)))
close(HeadFile)

##################################################################################
# GOEnrichment
GOEnrichment = fread(GOFile, h=T, quote="", sep="\t", select = c(1:4), stringsAsFactors=F)
GOEnrichment$PropObserved = GOEnrichment$ForegroundFrequency/ForegroundCount
GOEnrichment$PropExpected = GOEnrichment$BackgroundFrequency/BackgroundCount

# Hypergeometric test
GOEnrichment$HyperPval = apply(GOEnrichment, 1, function(x) phyper(as.numeric(x["ForegroundFrequency"])-1, as.numeric(x["BackgroundFrequency"]), 
                                                              BackgroundCount-as.numeric(x["BackgroundFrequency"]), ForegroundCount, lower.tail=FALSE))

# Binomial test
GOEnrichment$BinomPval = apply(GOEnrichment, 1, function(x) binom.test(as.numeric(x["ForegroundFrequency"]), ForegroundCount,
                                                                       as.numeric(x["BackgroundFrequency"])/BackgroundCount,
                                                                       alternative = "greater")$p.value)

# FDR                   
GOEnrichment$FDR = p.adjust(GOEnrichment$HyperPval, method = "fdr")
GOEnrichment$BinomFDR = p.adjust(GOEnrichment$BinomPval, method = "fdr")

##################################################################################
# Output
GOEnrichment <- GOEnrichment[order(GOEnrichment$FDR),]
write.table(GOEnrichment, file = paste0(pathResult, dist, "/FANTOM5.first1000.kidney.enhancers.hg38/GOEnrichment/", file),
            quote=FALSE, sep = '\t', dec = '.', row.names = F, col.names = T)

##################################################################################

