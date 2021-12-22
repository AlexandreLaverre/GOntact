path="/home/laverre/GOntact/results/"

sp = "human"
GOtype = "complete"
Trans =  "" #".Trans"
bait2bait = "" # ".bait2bait"
propagation = ".WithoutPropagation"

samples = c("kidney", "neuron", "liver", "heart")

GOEnrichment.list <- list()

for (sample in samples){
  print(sample)
  enh=paste("FANTOM5.first1000", sample, "enhancers.hg38", sep=".")
  
  GOFile = paste0(path, sp, "/", enh, "/", GOtype, ".GO.frequency", bait2bait, Trans, propagation, ".txt")
  
  ## Get the total Foreground & Background count in the header
  HeadFile <- file(GOFile,"r")
  ForegroundCount = as.numeric(sub(".*= ", "",  readLines(HeadFile,n=1)))
  BackgroundCount = as.numeric(sub(".*= ", "",  readLines(HeadFile,n=1)))
  close(HeadFile)
  

  GOEnrichment = read.table(GOFile, h=T, quote="", sep="\t")
  
  # Hypergeometric test
  GOEnrichment$Pval = apply(GOEnrichment, 1, function(x) phyper(as.numeric(x["ForegroundFrequency"])-1, as.numeric(x["BackgroundFrequency"]), 
                                                                BackgroundCount-as.numeric(x["BackgroundFrequency"]), ForegroundCount, lower.tail=FALSE))
  # FDR                   
  GOEnrichment$FDR = p.adjust(GOEnrichment$Pval, method = "fdr")
  GOEnrichment = GOEnrichment[which(GOEnrichment$FDR < 10e-3),]

  GOEnrichment.list[[sample]][["bio.process"]] = GOEnrichment[which(GOEnrichment$GONamespace == "biological_process"),]
  GOEnrichment.list[[sample]][["cell.compo"]] = GOEnrichment[which(GOEnrichment$GONamespace == "cellular_component"),]
  GOEnrichment.list[[sample]][["mol.funct"]] = GOEnrichment[which(GOEnrichment$GONamespace == "molecular_function"),]

}

#######################################################################################################
# Distance distribution between fragments

par(mfrow=c(2,2))
for (sample in samples){
  enh=paste("FANTOM5.first1000", sample, "enhancers.hg38", sep=".")
  
  ForegroundContact = read.table(paste0(path, sp, "/", enh, "/foreground.contacts.txt"), h=T, sep="\t")
  hist(ForegroundContact$Distance, breaks=50, main=sample, xlab="Distance")
  med = median(ForegroundContact$Distance)
  abline(v=med, col="red")
  mtext(paste0("med.dist=", med, "\n", "N=", nrow(ForegroundContact)), line = -1)
  
}