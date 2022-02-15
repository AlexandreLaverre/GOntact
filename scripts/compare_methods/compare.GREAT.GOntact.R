##################################################################################

pathGREAT="../../results/GO_enrichment_GREAT/"
pathGOntact="../../results/GO_enrichment_contacts/"
pathComparison="../../results/method_comparison/"
pathFigures="../../results/figures/"

##################################################################################

sp="mouse"
dataset="Maqbool2020_DN_mm10"
#dataset="VistaLimbEnhancerGenie_mm10"
background="ENCODE.Laverre2022"
space="biological_process"

maxFDR=0.05
minEnrichment=1.5

##################################################################################

great=read.table(paste(pathGREAT, sp,"/classical_upstream5kb_downstream1kb_extend1Mb/",dataset,"/enrichment_",space,"_Ensembl94_background",background,".txt",sep=""), h=T, stringsAsFactors=F, sep="\t", quote="\"")

great$Enrichment=great$Observed/great$Expected

rownames(great)=great$ID

gontact=read.table(paste(pathGOntact, sp, "/mindist0Kb_maxdist1Mb_extendOverlap5Kb/",dataset,"/GOEnrichment/biological_process.Background.",background,".txt", sep=""), h=T, stringsAsFactors=F, sep="\t", quote="\"")

rownames(gontact)=gontact$GOTerm

all.tested=intersect(rownames(gontact), rownames(great))
gontact=gontact[all.tested,]
great=great[all.tested,]

##################################################################################

signif.great=great$ID[which(great$FDR<maxFDR & great$Enrichment>=minEnrichment)]
signif.gontact=gontact$GOTerm[which(gontact$BinomFDR<maxFDR & gontact$Enrichment>=minEnrichment)]

onlygreat=setdiff(signif.great, signif.gontact)
onlygontact=setdiff(signif.gontact, signif.great)
common=intersect(signif.great, signif.gontact)

signif.any=unique(c(great$ID[which(great$FDR<maxFDR)], gontact$GOTerm[which(gontact$BinomFDR<maxFDR)]))

##################################################################################

pdf(file=paste(pathFigures, "GREAT_vs_GOntact_",sp,"_",dataset,".pdf",sep=""), width=6, height=6)

par(mar=c(4.1, 4.1, 2.1, 1.1))
plot(great[signif.any, "Enrichment"], gontact[signif.any, "Enrichment"], pch=20, xlab="", ylab="", main="", axes=F, col="gray70", xlim=c(0, 17), ylim=c(0, 17))
points(great[onlygreat, "Enrichment"], gontact[onlygreat, "Enrichment"], pch=20, col="red")
points(great[onlygontact, "Enrichment"], gontact[onlygontact, "Enrichment"], pch=20, col="blue")
points(great[common, "Enrichment"], gontact[common, "Enrichment"], pch=20, col="black")
axis(side=1, mgp=c(3, 0.75, 0))
axis(side=2, mgp=c(3, 0.75, 0))
mtext("enrichment, GREAT", side=1, line=2.5, cex=0.9)
mtext("enrichment, GOntact", side=2, line=2.5, cex=0.9)

legend("topleft", legend=c("only GREAT", "only GOntact", "both methods"), pch=20, col=c("red", "blue", "black"), inset=0.01, bg="white")

abline(0,1, lty=3, col="gray40")

dev.off()

##################################################################################

assoc.enh=read.table(paste(pathComparison, sp, "/",dataset,"/biological_process_GREAT.classical_upstream5kb_downstream1kb_extend1Mb_vs_GOntact.mindist0Kb_maxdist1Mb_extendOverlap5Kb.txt",sep=""), h=T, stringsAsFactors=F)

rownames(assoc.enh)=assoc.enh$GOTerm

assoc.enh$PropCommonMethod1=assoc.enh$NbInCommon/assoc.enh$NbEnhancersMethod1
assoc.enh$PropCommonMethod2=assoc.enh$NbInCommon/assoc.enh$NbEnhancersMethod2


signif.great.relaxed=great$ID[which(great$FDR<maxFDR)]
signif.gontact.relaxed=gontact$GOTerm[which(gontact$BinomFDR<maxFDR)]

common.relaxed=intersect(signif.great.relaxed, signif.gontact.relaxed)

##################################################################################


pdf(file=paste(pathFigures, "GREAT_vs_GOntact_nb_enhancers_",sp,"_",dataset,".pdf",sep=""), width=10, height=6)

par(mfrow=c(1,2))

ymax=max(c(assoc.enh[,"NbEnhancersMethod1"], assoc.enh[,"NbEnhancersMethod2"]))

plot(assoc.enh[,"NbEnhancersMethod1"], assoc.enh[,"NbEnhancersMethod2"], pch=20, xlab="", ylab="", axes=F, main="", ylim=c(0, ymax), xlim=c(0,ymax))

axis(side=1, mgp=c(3, 0.75, 0))
axis(side=2, mgp=c(3, 0.75, 0))
mtext("nb. associated enhancers, GREAT", side=1, line=2.5, cex=0.9)
mtext("nb. associated enhancers, GOntact", side=2, line=2.5, cex=0.9)

abline(0, 1, col="red", lty=3)

mtext("all GO terms", side=3)


ymax=max(c(assoc.enh[common.relaxed,"NbEnhancersMethod1"], assoc.enh[common.relaxed,"NbEnhancersMethod2"]))

plot(assoc.enh[common.relaxed,"NbEnhancersMethod1"], assoc.enh[common.relaxed,"NbEnhancersMethod2"], pch=20, xlab="", ylab="", axes=F, main="", ylim=c(0, ymax), xlim=c(0,ymax))


axis(side=1, mgp=c(3, 0.75, 0))
axis(side=2, mgp=c(3, 0.75, 0))
mtext("nb. associated enhancers, GREAT", side=1, line=2.5, cex=0.9)
mtext("nb. associated enhancers, GOntact", side=2, line=2.5, cex=0.9)

mtext("significant GO terms", side=3)

abline(0, 1, col="red", lty=3)

dev.off()

##################################################################################

pdf(file=paste(pathFigures, "GREAT_vs_GOntact_shared_enhancers_",sp,"_",dataset,".pdf",sep=""), width=3.5, height=6)
par(mar=c(4.1, 4.1, 2.1, 1.1))

boxplot(assoc.enh[common.relaxed, "PropCommonMethod1"], assoc.enh[common.relaxed, "PropCommonMethod2"], pch=20,  boxwex=0.5, names=rep("",2))
mtext(c("nb. shared/\nnb. GREAT", "nb.shared/\nnb. GOntact"), at=1:2, side=1, line=1.5, las=1, cex=0.9)
mtext("proportion shared enhancers, per GO term", side=2, line=2.5, cex=0.9)
dev.off()


##################################################################################
