##################################################################

pathGO="../../data/GeneOntology/"

##################################################################

gocat=read.table(paste(pathGO, "GOCategories.txt",sep=""), h=T, stringsAsFactors=F,sep="\t", quote="\"")

for(sp in c("human", "mouse")){
  for(space in c("biological_process", "molecular_function", "cellular_component")){
    this.gocat=gocat[which(gocat$GOSpace==space),]
    this.clusters=read.table(paste(pathGO, sp, ".GO.clusters.",space,".txt", sep=""), h=F, stringsAsFactors=F,sep="\t", quote="\"")

    ucl=unique(this.clusters$V1)
    desc=unlist(lapply(ucl, function(x) this.clusters$V3[which(this.clusters$V1==x)[1]]))

    this.res=data.frame("ID"=ucl, "GOSpace"=rep(space, length(ucl)), "Name"=desc, stringsAsFactors=F)
    
    this.gocat=this.gocat[which(!this.gocat$ID%in%this.clusters$V2),] ## remove clusters
    this.gocat=rbind(this.gocat, this.res, stringsAsFactors=F)

    write.table(this.gocat, paste(pathGO, "GOCategories.",sp,".",space,".txt",sep=""), row.names=F, col.names=T, sep="\t", quote=F)
  }

}

##################################################################
