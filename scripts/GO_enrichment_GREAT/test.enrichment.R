###########################################################################

pathResults="../../results/GO_enrichment_GREAT/"
pathEnsembl="../../data/ensembl_annotations/"
pathGO="../../data/GeneOntology/"

ensrelease=94
basal5=5000
basal3=1000
max.extend=1e+06

###########################################################################

sp="human"
dataset="FANTOM5.first1000.kidney.enhancers.hg38"

###########################################################################

go=read.table(paste(pathGO,"GOCategories.txt",sep=""),h=T, stringsAsFactors=F,sep="\t", quote="\"")
rownames(go)=go$ID

###########################################################################

chr.sizes=read.table(paste(pathEnsembl, sp,"/chr_sizes_Ensembl",ensrelease,".txt",sep=""),h=F, stringsAsFactors=F)
tot.size=sum(chr.sizes[,2])

###########################################################################

expected=read.table(paste(pathResults, sp, "/expected_values_Ensembl",ensrelease,"_basal5",basal5,"_basal3",basal3,"_max_extend",max.extend,".txt",sep=""), h=T, stringsAsFactors=F)

observed=read.table(paste(pathResults, sp, "/",dataset,"/observed_values_Ensembl",ensrelease,"_basal5",basal5,"_basal3",basal3,"_max_extend",max.extend,".txt",sep=""), h=T, stringsAsFactors=F)

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
    p=expected$NbBases[which(expected$ID==cat)]/tot.size

    pval=binom.test(x=x, n=n, p=p)$p.value

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

  write.table(res, file=paste(pathResults, sp,"/",dataset,"/enrichment_",space,"_Ensembl",ensrelease,"_basal5",basal5,"_basal3",basal3,"_max_extend",max.extend,".txt",sep=""), row.names=F, col.names=T, quote=F)
}

###########################################################################
