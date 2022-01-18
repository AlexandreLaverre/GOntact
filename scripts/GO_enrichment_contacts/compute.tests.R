##################################################################################
#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)

sp= args[1]
enhancer= args[2]
dist = args[3]
file = args[4]

message("### Perform GOEnrichment test ###")
pathResult = "/beegfs/data/necsulea/GOntact/results/GO_enrichment_contacts/"
GOFile = paste(pathResult, sp, dist, enhancer, "GOFrequency", file, sep="/")

##################################################################################
## Get the total Foreground & Background count in the header
HeadFile <- file(GOFile,"r")
ForegroundCount = as.numeric(sub(".*= ", "",  readLines(HeadFile,n=1)))
BackgroundCount = as.numeric(sub(".*= ", "",  readLines(HeadFile,n=1)))
close(HeadFile)

##################################################################################
# GOEnrichment
GOEnrichment = read.table(GOFile, h=T, quote="", sep="\t", select = c(1:4), stringsAsFactors=F)
GOEnrichment$PropObserved = GOEnrichment$ForegroundFrequency/ForegroundCount
GOEnrichment$PropExpected = GOEnrichment$BackgroundFrequency/BackgroundCount
GOEnrichment$Enrichment=GOEnrichment$PropObserved/GOEnrichment$PropExpected

# Hypergeometric test
GOEnrichment$HyperPval = apply(GOEnrichment, 1, function(x) phyper(as.numeric(x["ForegroundFrequency"])-1, as.numeric(x["BackgroundFrequency"]), 
                                                              BackgroundCount-as.numeric(x["BackgroundFrequency"]), ForegroundCount, lower.tail=FALSE))

# Binomial test
GOEnrichment$BinomPval = apply(GOEnrichment, 1, function(x) binom.test(as.numeric(x["ForegroundFrequency"]), ForegroundCount,
                                                                       as.numeric(x["BackgroundFrequency"])/BackgroundCount,
                                                                       alternative = "greater")$p.value)

# FDR                   
GOEnrichment$HyperFDR = p.adjust(GOEnrichment$HyperPval, method = "fdr")
GOEnrichment$BinomFDR = p.adjust(GOEnrichment$BinomPval, method = "fdr")

##################################################################################
# Output
GOEnrichment <- GOEnrichment[order(GOEnrichment$HyperFDR),]
write.table(GOEnrichment, file=paste(pathResult, sp, dist, enhancer, "GOEnrichment", file, sep="/"),
            quote=FALSE, sep = '\t', dec = '.', row.names = F, col.names = T)

##################################################################################

