########################################################################

pathEnhancers="../../data/enhancers/"
pathFigures="../../data_for_publication/figures/"

options(stringsAsFactors=F)

########################################################################

for(sp in c("human", "mouse")){
  encode=read.table(paste(pathEnhancers, sp, "/ENCODE.Laverre2022.bed",sep=""), h=F, stringsAsFactors=F)

  save(encode, file=paste(pathFigures, "RData/data.ENCODE.",sp,".RData",sep=""))
}

########################################################################
