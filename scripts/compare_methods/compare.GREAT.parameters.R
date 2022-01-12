##################################################################################

pathResults="../../results/GO_enrichment_GREAT/"
pathFigures="../../results/figures/"

##################################################################################

sp="human"
dataset="FANTOM5.first1000.kidney.enhancers.hg38"
space="biological_process"
background="regulatory_regions"

method1="classical_upstream5kb_downstream1kb_extend1Mb"
method2="classical_mindist25kb_extend1Mb"

maxFDR=0.01

##################################################################################

## GO enrichment tests

res1=read.table(paste(pathResults, sp,"/",method1,"/",dataset,"/enrichment_regulatory_regions_",space,"_Ensembl94.txt",sep=""), h=T, stringsAsFactors=F, sep="\t", quote="\"")

rownames(res1)=res1$ID

res2=read.table(paste(pathResults, sp,"/",method2,"/",dataset,"/enrichment_regulatory_regions_",space,"_Ensembl94.txt",sep=""), h=T, stringsAsFactors=F, sep="\t", quote="\"")

rownames(res2)=res2$ID

##################################################################################

nbsignif1=length(which(res1$FDR<maxFDR))
nbsignif2=length(which(res2$FDR<maxFDR))

common=intersect(res1$ID[which(res1$FDR<maxFDR)], res2$ID[which(res2$FDR<maxFDR)])
any=union(res1$ID[which(res1$FDR<maxFDR)], res2$ID[which(res2$FDR<maxFDR)])

only1=setdiff(res1$ID[which(res1$FDR<maxFDR)], res2$ID[which(res2$FDR<maxFDR)])
only2=setdiff(res2$ID[which(res2$FDR<maxFDR)], res1$ID[which(res1$FDR<maxFDR)])

##################################################################################

## compare enrichment

res1$Enrichment=res1$Observed/res1$Expected
res2$Enrichment=res2$Observed/res2$Expected

lim=range(c(res1[any, "Enrichment"],res2[any, "Enrichment"]))

pdf(paste(pathFigures, "enrichment_",method1,"_vs_",method2,".pdf",sep=""))

plot(res1[any, "Enrichment"], res2[any, "Enrichment"], pch=20, xlab="", ylab="", axes=F, main="", xlim=lim, ylim=lim)
points(res1[common, "Enrichment"], res2[common, "Enrichment"], pch=20, col="blue")
points(res1[only1, "Enrichment"], res2[only1, "Enrichment"], pch=20, col="red")
points(res1[only2, "Enrichment"], res2[only2, "Enrichment"], pch=20, col="orange")

abline(0,1)

axis(side=1, mgp=c(3, 0.5, 0))
mtext(method1, side=1, line=2)

axis(side=2, mgp=c(3, 0.75, 0))
mtext(method2, side=2, line=2)

dev.off()

##################################################################################
