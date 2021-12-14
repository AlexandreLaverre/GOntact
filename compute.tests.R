sp = "human"
enhancers = "FANTOM5.brain.enhancers.hg38"
GOtype = "complete"

## Get the total Foreground & Background count in the header
HeadFile <- file(paste0("/home/laverre/GOntact/results/", sp, "/", enhancers, "/", GOtype, ".GO.frequency.txt"),"r")
ForegroundCount = as.numeric(sub(".*= ", "",  readLines(HeadFile,n=1)))
BackgroundCount = as.numeric(sub(".*= ", "",  readLines(HeadFile,n=1)))
close(HeadFile)

GOEnrichment = read.table(paste0("/home/laverre/GOntact/results/", sp, "/", enhancers, "/", GOtype, ".GO.frequency.txt"), h=T, quote="", sep="\t")

# Hypergeometric test
GOEnrichment$Pval = apply(GOEnrichment, 1, function(x) phyper(as.numeric(x["ForegroundFrequency"]), as.numeric(x["BackgroundFrequency"]), 
                                                              BackgroundCount-as.numeric(x["BackgroundFrequency"]), ForegroundCount, lower.tail=FALSE))
# FDR                   
GOEnrichment$FDR = p.adjust(GOEnrichment$Pval, method = "fdr")


