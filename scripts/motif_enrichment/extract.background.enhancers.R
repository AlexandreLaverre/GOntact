#####################################################################################

pathResults="../../results/motif_enrichment/"
pathEnhancers="../../data/enhancers/"

options(stringsAsFactors=F)

#####################################################################################

n=50000

for(sp in c("human", "mouse")){
  enhancers=read.table(paste(pathEnhancers, sp, "/ENCODE.Laverre2022.bed",sep=""), h=F, stringsAsFactors=F)

  selected=sample(1:nrow(enhancers), size=n, re=F)

  write.table(enhancers[selected,], file=paste(pathResults,sp, "/ENCODE.Laverre2022.random50K.bed",sep=""), row.names=F, col.names=F, sep="\t", quote=F)
}

#####################################################################################
