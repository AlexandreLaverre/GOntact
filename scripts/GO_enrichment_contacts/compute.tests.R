sp = "human"
enhancers = "FANTOM5.first150.kidney.enhancers.hg38" #"FANTOM5.first150.neuron.enhancers.hg38" 
GOtype = "complete"
Trans =  "" #".Trans"
bait = "" # ".bait2bait"

GOEnrichment.list <- list()

for (propagation in c("", "before.propagation")){
  for (Baited in c("", ".BaitedEnh")){
    
    GOFile = paste0("/home/laverre/GOntact/results/", sp, "/", propagation, "/", enhancers, "/", GOtype, ".GO.frequency", Baited, Trans, bait, ".txt")
  
    ## Get the total Foreground & Background count in the header
    HeadFile <- file(GOFile,"r")
    ForegroundCount = as.numeric(sub(".*= ", "",  readLines(HeadFile,n=1)))
    BackgroundCount = as.numeric(sub(".*= ", "",  readLines(HeadFile,n=1)))
    close(HeadFile)
    
    print(GOFile)
    GOEnrichment = read.table(GOFile, h=T, quote="", sep="\t")
    
    # Hypergeometric test
    GOEnrichment$Pval = apply(GOEnrichment, 1, function(x) phyper(as.numeric(x["ForegroundFrequency"]), as.numeric(x["BackgroundFrequency"]), 
                                                                  BackgroundCount-as.numeric(x["BackgroundFrequency"]), ForegroundCount, lower.tail=FALSE))
    # FDR                   
    GOEnrichment$FDR = p.adjust(GOEnrichment$Pval, method = "fdr")
    GOEnrichment = GOEnrichment[which(GOEnrichment$FDR < 10e-3),]
    
    propagation.name = ifelse(propagation == "", "after.propagation", "before.propagation")
    Baited.name = ifelse(Baited == "", "default", "baited.enh")
    type = paste(propagation.name, Baited.name, sep=".")
    
    GOEnrichment.list[[paste0(type, ".bio")]] = GOEnrichment[which(GOEnrichment$GONamespace == "biological_process"),]
    GOEnrichment.list[[paste0(type,".cell")]] = GOEnrichment[which(GOEnrichment$GONamespace == "cellular_component"),]
    GOEnrichment.list[[paste0(type,".mol")]] = GOEnrichment[which(GOEnrichment$GONamespace == "molecular_function"),]
    
    
  }
}

#######################################################################################################

View(GOEnrichment.list[["before.propagation.default.bio"]])
View(GOEnrichment.list[["after.propagation.default.bio"]])

View(GOEnrichment.list[["before.propagation.baited.enh.bio"]])
View(GOEnrichment.list[["after.propagation.baited.enh.bio"]])


#######################################################################################################