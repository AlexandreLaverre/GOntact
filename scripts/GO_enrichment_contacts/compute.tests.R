sp = "human"
enhancers = "FANTOM5.first150.Tcell.hg38" #"FANTOM5.first150.neuron.enhancers.hg38" 
GOtype = "complete"
Baited = ".BaitedEnh"
Trans =  "" #".Trans"
bait = "" # ".bait2bait"

## Get the total Foreground & Background count in the header
HeadFile <- file(paste0("/home/laverre/GOntact/results/", sp, "/", enhancers, "/", GOtype, ".GO.frequency", Baited, Trans, bait, ".txt"),"r")
ForegroundCount = as.numeric(sub(".*= ", "",  readLines(HeadFile,n=1)))
BackgroundCount = as.numeric(sub(".*= ", "",  readLines(HeadFile,n=1)))
close(HeadFile)

GOEnrichment = read.table(paste0("/home/laverre/GOntact/results/", sp, "/", enhancers, "/", GOtype, ".GO.frequency", Baited, Trans, bait, ".txt"), h=T, quote="", sep="\t")

# Hypergeometric test
GOEnrichment$Pval = apply(GOEnrichment, 1, function(x) phyper(as.numeric(x["ForegroundFrequency"]), as.numeric(x["BackgroundFrequency"]), 
                                                              BackgroundCount-as.numeric(x["BackgroundFrequency"]), ForegroundCount, lower.tail=FALSE))
# FDR                   
GOEnrichment$FDR = p.adjust(GOEnrichment$Pval, method = "fdr")


GOEnrichment.default = GOEnrichment
