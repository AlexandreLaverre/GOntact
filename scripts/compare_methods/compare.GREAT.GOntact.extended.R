##################################################################################

pathGREAT="../../results/GO_enrichment_GREAT/"
pathGOntact="../../results/GO_enrichment_contacts/"
pathComparison="../../results/method_comparison/"
pathFigures="../../results/figures/"

##################################################################################

sp="mouse"
samples=c("Maqbool2020_DP_mm10", "Maqbool2020_DN_mm10","Maqbool2020_Th2_mm10")
synonyms=c("DP T cells", "DN T cells", "Th2 T cells")
background="ENCODE.Laverre2022"
space="biological_process"

methods.great=c("classical_upstream5kb_downstream1kb_extend1Mb", "classical_mindist25kb_extend1Mb")
methods.gontact=c("mindist0Kb_maxdist1Mb_extendOverlap5Kb", "mindist25Kb_maxdist1Mb_extendOverlap5Kb")

nb.methods=2

maxFDR=0.05
minEnrichment=1.5

##################################################################################

mat.signif.great=matrix(rep(NA, length(samples)*nb.methods), nrow=nb.methods)
rownames(mat.signif.great)=c("maxDist0Kb_minDist1Mb", "maxDist25Kb_minDist1Mb")
colnames(mat.signif.great)=samples

mat.signif.gontact=mat.signif.great
mat.signif.common=mat.signif.great

##################################################################################

for(sample in samples){
  for(i in 1:nb.methods){
    this.great=methods.great[i]
    this.gontact=methods.gontact[i]
    
    great=read.table(paste(pathGREAT, sp,"/",this.great,"/",sample,"/enrichment_",space,"_Ensembl94_background",background,".txt",sep=""), h=T, stringsAsFactors=F, sep="\t", quote="\"")

    great$Enrichment=great$Observed/great$Expected
    
    rownames(great)=great$ID
    
    gontact=read.table(paste(pathGOntact, sp, "/",this.gontact,"/",sample,"/GOEnrichment/biological_process.Background.",background,".txt", sep=""), h=T, stringsAsFactors=F, sep="\t", quote="\"")
    
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

    mat.signif.great[i,sample]=length(signif.great)
    mat.signif.gontact[i,sample]=length(signif.gontact)
    mat.signif.common[i,sample]=length(common)
    
   ##################################################################################
  }
}

## assoc.enh=read.table(paste(pathComparison, sp, "/",dataset,"/biological_process_GREAT.classical_upstream5kb_downstream1kb_extend1Mb_vs_GOntact.mindist0Kb_maxdist1Mb_extendOverlap5Kb.txt",sep=""), h=T, stringsAsFactors=F)

## rownames(assoc.enh)=assoc.enh$GOTerm

## assoc.enh$PropCommonMethod1=assoc.enh$NbInCommon/assoc.enh$NbEnhancersMethod1
## assoc.enh$PropCommonMethod2=assoc.enh$NbInCommon/assoc.enh$NbEnhancersMethod2

##################################################################################

method.xpos=c(0, 0.5) ## two methods
tool.xpos=c(0, 1.5, 3) ## great, gontact, common
names(tool.xpos)=c("great","gontact","common")

sample.xpos=c(0, 6, 12)
names(sample.xpos)=samples

width=0.15

ymax=max(c(as.numeric(mat.signif.great), as.numeric(mat.signif.gontact), as.numeric(mat.signif.common)))

colors=c("gray20", "gray60")

pdf(file=paste(pathFigures, "GREAT_vs_GOntact_",sp,"_extended.pdf",sep=""), width=6, height=6)

par(mar=c(6.1, 4.1, 2.1, 1.1))
plot(1, type="n", xlim=c(-1, 17), ylim=c(0,ymax+10), ylab="",xlab="", axes=F)
axis(side=2)
mtext("number of significant GO terms (FDR < 0.05, enrichment > 1.5)", side=2, line=2.5, cex=0.95)

axis(side=1, at=kronecker(as.numeric(sample.xpos),c(0.25, 1.75, 3.25),function(x,y) {x+y}), labels=rep("", 9))

for(sample in samples){
for(tool in c("great", "gontact", "common")){
  mat=get(paste("mat.signif.",tool,sep=""))
  for(i in 1:2){
    this.x=tool.xpos[tool]+sample.xpos[sample]+method.xpos[i]

      this.y=mat[i,sample]

      rect(this.x-width,0, this.x+width,this.y,col=colors[i],bg=colors[i])

      print(paste(tool, sample, i, this.x))
      
    }
  }
}

legend("topleft", legend=c("0kb - 1Mb", "25kb - 1Mb"), fill=c("gray20", "gray60"), inset=c(0.01, 0.1))
mtext(rep(c("GREAT","GOntact", "common"), 3), side=1, at=kronecker(as.numeric(sample.xpos),c(0.25, 1.75, 3.25),function(x,y) {x+y}), line=0.75, las=2, cex=0.95)

mtext(synonyms[1], side=1, at=sample.xpos[1]+1.5, line=5)
mtext(synonyms[2], side=1, at=sample.xpos[2]+1.5, line=5)
mtext(synonyms[3], side=1, at=sample.xpos[3]+1.5, line=5)

dev.off()

##################################################################################
