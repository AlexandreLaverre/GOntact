###########################################################################

pathResults="../../results/GO_enrichment_GREAT/"
pathGO="../../data/GeneOntology/"

ensrelease=94
basal5=5000
basal3=1000
max.extend=1e+06
outdir="classical_upstream5kb_downstream1kb_extend1Mb"

###########################################################################

sp="human"
dataset="FANTOM5.first1000.kidney.enhancers.hg38"

###########################################################################

go=read.table(paste(pathGO,"GOCategories.txt",sep=""),h=T, stringsAsFactors=F,sep="\t", quote="\"")
rownames(go)=go$ID

###########################################################################

info.chr=readLines(paste(pathResults, sp, "/",outdir,"/expected_values_Ensembl",ensrelease,".txt",sep=""))[1:2]
  
totsize=as.numeric(unlist(strsplit(info.chr[1], split="\t"))[2])
nbN=as.numeric(unlist(strsplit(info.chr[2], split="\t"))[2])

size.nonN=totsize-nbN

###########################################################################

expected=read.table(paste(pathResults, sp, "/",outdir,"/expected_values_Ensembl",ensrelease,".txt",sep=""), h=T, stringsAsFactors=F)

observed=read.table(paste(pathResults, sp, "/",outdir,"/",dataset,"/observed_values_Ensembl",ensrelease,".txt",sep=""), h=T, stringsAsFactors=F)

###########################################################################

for(space in unique(expected$GOSpace)){
  print(space)
  
  tested.cat=c()
  tested.observed=c()
  tested.expected=c()
  tested.pval=c()
  
  for(cat in expected$ID[which(expected$GOSpace==space)]){
    x=observed$NbAssociatedElements[which(observed$ID==cat)]
    n=observed$NbTotalElements[which(observed$ID==cat)]
    p=expected$NbNonNBases[which(expected$ID==cat)]/size.nonN

    pval=1-pbinom(q=x/n, size=n, p=p)

    tested.cat=c(tested.cat, cat)
    tested.observed=c(tested.observed, x/n)
    tested.expected=c(tested.expected, p)
    tested.pval=c(tested.pval, pval)
  }

  res=data.frame("ID"=tested.cat, "Observed"=tested.observed,"Expected"=tested.expected, "PValue"=tested.pval, stringsAsFactors=F)

  res$Name=go[res$ID,"Name"]

  res$FDR=p.adjust(res$PValue)

  res=res[order(res$FDR),]

  res=res[,c("ID","Name", "Observed","Expected", "PValue", "FDR")]
  res$Expected=round(res$Expected, digits=2)
  res$Observed=round(res$Observed, digits=2)

  write.table(res, file=paste(pathResults, sp,"/",outdir, "/",dataset,"/enrichment_",space,"_Ensembl",ensrelease,".txt",sep=""), row.names=F, col.names=T, quote=F)
}

###########################################################################
