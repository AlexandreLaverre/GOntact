##################################################################################
#!/usr/bin/env Rscript

path = "/beegfs/data/necsulea/GOntact/" 
pathResults= paste0(path, "/results/GO_enrichment_contacts/")
pathFigures= paste0(path, "/results/GO_enrichment_contacts/figures/")
args = commandArgs(trailingOnly=TRUE)

##################################################################################
#args = commandArgs(trailingOnly=TRUE)

sp="human"
dataset="FANTOM5.first1000.kidney.enhancers.hg38"
space="biological_process"

method1=args[1] # i.e "mindist0Kb_maxdist2Mb_extendOverlap5Kb"
method2=args[2] # i.e : "mindist0Kb_maxdist1Mb_extendOverlap5Kb"
background=args[3] # i.e : FANTOM5.Laverre2022 or FANTOM5.selection

maxFDR=0.05

##################################################################################
## GO enrichment tests

res1=read.table(paste0(pathResults, sp,"/",method1,"/",dataset,"/GOEnrichment/uniqueGO.",space,".Background.EnhancersContacts",
                      background, ".EnhancerScale.EnhancerCountOnce.txt"), h=T, stringsAsFactors=F, sep="\t", quote="\"")

rownames(res1)=res1$GOTerm

res2=read.table(paste0(pathResults, sp,"/",method2,"/",dataset,"/GOEnrichment/uniqueGO.",space,".Background.EnhancersContacts",
                       background, ".EnhancerScale.EnhancerCountOnce.txt"), h=T, stringsAsFactors=F, sep="\t", quote="\"")
rownames(res2)=res2$GOTerm

##################################################################################

nbsignif1=length(which(res1$FDR<maxFDR))
nbsignif2=length(which(res2$FDR<maxFDR))

common=intersect(res1$GOTerm[which(res1$FDR<maxFDR)], res2$GOTerm[which(res2$FDR<maxFDR)])
any=union(res1$GOTerm[which(res1$FDR<maxFDR)], res2$GOTerm[which(res2$FDR<maxFDR)])

only1=setdiff(res1$GOTerm[which(res1$FDR<maxFDR)], res2$GOTerm[which(res2$FDR<maxFDR)])
only2=setdiff(res2$GOTerm[which(res2$FDR<maxFDR)], res1$GOTerm[which(res1$FDR<maxFDR)])

##################################################################################
## compare enrichment
res1$Enrichment=res1$PropObserved/res1$PropExpected
res2$Enrichment=res2$PropObserved/res2$PropExpected

lim=range(c(res1[any, "Enrichment"],res2[any, "Enrichment"]))

pdf(paste(pathFigures, "GOntact_enrichment_",method1,"_vs_",method2,"_backgound_",background, ".pdf",sep=""))

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
