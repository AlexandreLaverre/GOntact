###########################################################################

pathResults="../../results/GO_enrichment_GREAT/"
pathGO="../../data/GeneOntology/"

ensrelease=94

###########################################################################

args <- commandArgs(trailingOnly = TRUE)

sp=args[1]
dataset=args[2]
method=args[3]

print(paste("sp=",sp)) ## human or mouse
print(paste("dataset=",dataset)) ## e.g. "FANTOM5.first1000.kidney.enhancers.hg38"
print(paste("method=",method)) ## e.g. "classical_upstream5kb_downstream1kb_extend1Mb"

###########################################################################

go=read.table(paste(pathGO,"GOCategories.txt",sep=""),h=T, stringsAsFactors=F,sep="\t", quote="\"")
rownames(go)=go$ID

###########################################################################

info.chr=readLines(paste(pathResults, sp, "/",method,"/expected_values_Ensembl",ensrelease,".txt",sep=""))[1:4]

## total chromosome size
totsize=as.numeric(unlist(strsplit(info.chr[1], split="\t"))[2])
nbN=as.numeric(unlist(strsplit(info.chr[2], split="\t"))[2])

genomesize.nonN=totsize-nbN

print(paste("genome size without Ns", genomesize.nonN))

## ungapped region size

regionsize.nonN=as.numeric(unlist(strsplit(info.chr[4], split="\t"))[2])

print(paste("region size without Ns", regionsize.nonN))

###########################################################################

expected=read.table(paste(pathResults, sp, "/",method,"/expected_values_Ensembl",ensrelease,".txt",sep=""), h=T, stringsAsFactors=F)

observed=read.table(paste(pathResults, sp, "/",method,"/",dataset,"/observed_values_Ensembl",ensrelease,".txt",sep=""), h=T, stringsAsFactors=F)

info=readLines(paste(pathResults, sp, "/",method,"/",dataset,"/observed_values_Ensembl",ensrelease,".txt",sep=""))[1:2]
ntot=as.numeric(unlist(strsplit(info[1],split="\t"))[2]) ## total number of elements
nreg=as.numeric(unlist(strsplit(info[2],split="\t"))[2]) ## total number of elements found in regulatory regions

print(paste(ntot,"elements in total"))
print(paste(nreg,"elements that overlap with regulatory regions"))

###########################################################################

for(space in unique(expected$GOSpace)){
  print(space)
  
  genome.tested.cat=c()
  genome.tested.observed=c()
  genome.tested.expected=c()
  genome.tested.pval=c()

  region.tested.cat=c()
  region.tested.observed=c()
  region.tested.expected=c()
  region.tested.pval=c()
    
  for(cat in intersect(expected$ID[which(expected$GOSpace==space)], go$ID)){
    p.genome=expected$NbNonNBases[which(expected$ID==cat)]/genomesize.nonN
    p.region=expected$NbNonNBases[which(expected$ID==cat)]/regionsize.nonN

    x=observed$NbAssociatedElements[which(observed$ID==cat)]
       
    pval.genome=binom.test(x, ntot, p=p.genome, alternative="greater")$p.value
    pval.region=binom.test(x, nreg, p=p.region, alternative="greater")$p.value

    genome.tested.cat=c(genome.tested.cat, cat)
    genome.tested.observed=c(genome.tested.observed, x/ntot)
    genome.tested.expected=c(genome.tested.expected, p.genome)
    genome.tested.pval=c(genome.tested.pval, pval.genome)

    region.tested.cat=c(region.tested.cat, cat)
    region.tested.observed=c(region.tested.observed, x/nreg)
    region.tested.expected=c(region.tested.expected, p.region)
    region.tested.pval=c(region.tested.pval, pval.region)
  }

  res.genome=data.frame("ID"=genome.tested.cat, "Observed"=genome.tested.observed,"Expected"=genome.tested.expected, "PValue"=genome.tested.pval, stringsAsFactors=F)
  res.genome$Name=go[res.genome$ID,"Name"]
  res.genome$FDR=p.adjust(res.genome$PValue)
  res.genome=res.genome[order(res.genome$FDR),]

  res.genome=res.genome[,c("ID","Name", "Observed","Expected", "PValue", "FDR")]

  write.table(res.genome, file=paste(pathResults, sp,"/",method, "/",dataset,"/enrichment_whole_genome_",space,"_Ensembl",ensrelease,".txt",sep=""), row.names=F, col.names=T, quote=F)


  ## only regulatory regions
  res.region=data.frame("ID"=region.tested.cat, "Observed"=region.tested.observed,"Expected"=region.tested.expected, "PValue"=region.tested.pval, stringsAsFactors=F)
  res.region$Name=go[res.region$ID,"Name"]
  res.region$FDR=p.adjust(res.region$PValue)
  res.region=res.region[order(res.region$FDR),]

  res.region=res.region[,c("ID","Name", "Observed","Expected", "PValue", "FDR")]

  write.table(res.region, file=paste(pathResults, sp,"/",method, "/",dataset,"/enrichment_regulatory_regions_",space,"_Ensembl",ensrelease,".txt",sep=""), row.names=F, col.names=T, quote=F)
}

###########################################################################
