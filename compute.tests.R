## Get the total Foreground & Background count in the header
HeadFile <- file("/home/laverre/GOntact/results/human/unique.GO.frequency.txt","r")
ForegroundCount = as.numeric(sub(".*= ", "",  readLines(HeadFile,n=1)))
BackgroundCount = as.numeric(sub(".*= ", "",  readLines(HeadFile,n=1)))
close(HeadFile)

GOEnrichment = read.table("/home/laverre/GOntact/results/human/unique.GO.frequency.txt", h=T)

# Hypergeometric test
GOEnrichment$Pval = apply(GOEnrichment, 1, function(x) phyper(as.numeric(x["ForegroundFrequency"]), as.numeric(x["BackgroundFrequency"]), 
                                                              BackgroundCount-as.numeric(x["BackgroundFrequency"]), ForegroundCount, lower.tail=FALSE))
# FDR                   
GOEnrichment$FDR = p.adjust(GOEnrichment$Pval, method = "fdr")


