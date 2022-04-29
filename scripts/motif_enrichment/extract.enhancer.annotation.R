#####################################################################################

pathResults="../../results/"

sample="Vista_heart_vs_ENCODE"

options(stringsAsFactors=F)

#####################################################################################

categories=c("GO:0060173")

#####################################################################################

for(sp in c("human", "mouse")){

  ## great
  great=read.table(paste(pathResults, sp, "/",sample, "/biological_process/GREAT_upstream5kb_downstream1kb_extend1Mb/GOntact_element_GO_association_background.txt", sep=""), h=F, stringsAsFactors=F)

  ## contacts
  gontact=read.table(paste(pathResults, sp, "/", sample, "/biological_process/contacts_mindist0kb_maxdist1Mb/GOntact_element_GO_association_background.txt",sep=""), h=F, stringsAsFactors=F)

  for(cat in categories){
    outcat=paste(unlist(strsplit(cat, split=":")),collapse="")

    system(paste("mkdir -p ",pathResults, "motif_enrichment/",sp,"/",outcat,"/GREAT_only/",sep=""))
    system(paste("mkdir -p ",pathResults, "motif_enrichment/",sp,"/",outcat,"/GOntact_only/",sep=""))
    system(paste("mkdir -p ",pathResults, "motif_enrichment/",sp,"/",outcat,"/shared/",sep=""))

    enh.great=great$V2[which(great$V1==cat)]
    enh.gontact=gontact$V2[which(gontact$V1==cat)]

    enh.great=unlist(strsplit(enh.great, split=","))
    enh.gontact=unlist(strsplit(enh.gontact, split=","))

    only.great=setdiff(enh.great, enh.gontact)
    only.gontact=setdiff(enh.gontact, enh.great)
    both=intersect(enh.gontact, enh.great)

    ## write output for only great
    chr.great=unlist(lapply(only.great, function(x) unlist(strsplit(x, split=":"))[1]))
    coords.great=unlist(lapply(only.great, function(x) unlist(strsplit(x, split=":"))[2]))
    start.great=unlist(lapply(coords.great, function(x) unlist(strsplit(x, split="-"))[1]))
    end.great=unlist(lapply(coords.great, function(x) unlist(strsplit(x, split="-"))[2]))

    res.great=data.frame(chr.great, start.great, end.great, only.great, rep(".", length(only.great)),rep("+", length(only.great)))

    write.table(res.great, paste(pathResults, "motif_enrichment/",sp,"/",outcat,"/GREAT_only/associated_enhancers.bed",sep=""), row.names=F, col.names=F, sep="\t", quote=F)

    ## write output for only gontact

    chr.gontact=unlist(lapply(only.gontact, function(x) unlist(strsplit(x, split=":"))[1]))
    coords.gontact=unlist(lapply(only.gontact, function(x) unlist(strsplit(x, split=":"))[2]))
    start.gontact=unlist(lapply(coords.gontact, function(x) unlist(strsplit(x, split="-"))[1]))
    end.gontact=unlist(lapply(coords.gontact, function(x) unlist(strsplit(x, split="-"))[2]))

    res.gontact=data.frame(chr.gontact, start.gontact, end.gontact, only.gontact, rep(".", length(only.gontact)),rep("+", length(only.gontact)))

    write.table(res.gontact, paste(pathResults, "motif_enrichment/",sp,"/",outcat,"/GOntact_only/associated_enhancers.bed",sep=""), row.names=F, col.names=F, sep="\t", quote=F)

    ## write output for shared

    chr.shared=unlist(lapply(both, function(x) unlist(strsplit(x, split=":"))[1]))
    coords.shared=unlist(lapply(both, function(x) unlist(strsplit(x, split=":"))[2]))
    start.shared=unlist(lapply(coords.shared, function(x) unlist(strsplit(x, split="-"))[1]))
    end.shared=unlist(lapply(coords.shared, function(x) unlist(strsplit(x, split="-"))[2]))

    res.shared=data.frame(chr.shared, start.shared, end.shared, both, rep(".", length(both)),rep("+", length(both)))

    write.table(res.shared, paste(pathResults, "motif_enrichment/",sp,"/",outcat,"/shared/associated_enhancers.bed",sep=""), row.names=F, col.names=F, sep="\t", quote=F)

  }

}

#####################################################################################
