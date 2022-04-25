#####################################################################################

pathFigures="../../data_for_publication/figures/"

options(stringsAsFactors=F)

#####################################################################################

if(load){
  load(paste(pathFigures, "RData/data.enrichment.results.RData",sep=""))

  species=c("human", "mouse")
  samples=c("Vista_heart_vs_ENCODE", "Vista_limb_vs_ENCODE")
  domain="biological_process"

  minFDR=1e-10
  maxFDR=0.01

  load=FALSE
}

#####################################################################################

if(prepare){

  results=list()

  for(i in 1:2){
    sp=species[i]
    sample=samples[i]

    ## enrichment results
    great1Mb=enrichment.results[[sp]][[sample]][[domain]][["GREAT_upstream5kb_downstream1kb_extend1Mb"]]
    rownames(great1Mb)=great1Mb$GOID
    great1Mb$Enrichment=great1Mb$Observed/great1Mb$Expected

    contacts1Mb=enrichment.results[[sp]][[sample]][[domain]][["contacts_mindist0kb_maxdist1Mb"]]
    rownames(contacts1Mb)=contacts1Mb$GOID
    contacts1Mb$Enrichment=contacts1Mb$Observed/contacts1Mb$Expected

    hybrid1Mb=enrichment.results[[sp]][[sample]][[domain]][["hybrid_mindist25kb_maxdist1Mb"]]
    rownames(hybrid1Mb)=hybrid1Mb$GOID
    hybrid1Mb$Enrichment=hybrid1Mb$Observed/hybrid1Mb$Expected

    ## replace very small FDRs

    great1Mb$FDR[which(great1Mb$FDR<minFDR)]=minFDR
    contacts1Mb$FDR[which(contacts1Mb$FDR<minFDR)]=minFDR
    hybrid1Mb$FDR[which(hybrid1Mb$FDR<minFDR)]=minFDR

    ## signif categories
    signif.great1Mb=great1Mb$GOID[which(great1Mb$FDR<maxFDR)]
    signif.contacts1Mb=contacts1Mb$GOID[which(contacts1Mb$FDR<maxFDR)]
    signif.hybrid1Mb=hybrid1Mb$GOID[which(hybrid1Mb$FDR<maxFDR)]
    signif.any1Mb=unique(c(signif.great1Mb, signif.contacts1Mb))

    ## direct comparison great vs contacts, 1Mb
    common1Mb=intersect(great1Mb$GOID, contacts1Mb$GOID)


    ## smaller names

    great1Mb$GOName[which(great1Mb$GOName=="positive regulation of epithelial cell proliferation")]="positive reg. of epithelial cell proliferation"
    great1Mb$GOName[which(great1Mb$GOName=="regulation of cellular macromolecule biosynthetic process")]="reg. of cell. macromolecule biosynthetic proc."
    contacts1Mb$GOName[which(contacts1Mb$GOName=="positive regulation of transcription by RNA polymerase II")]="positive reg. of transcription by RNA pol II"

    ## save results

    results[[sp]]=list("great1Mb"=great1Mb, "contacts1Mb"=contacts1Mb, "signif.great1Mb"=signif.great1Mb, "signif.contacts1Mb"=signif.contacts1Mb, "signif.any1Mb"=signif.any1Mb, "common1Mb")

  }

  prepare=FALSE
}

#####################################################################################

## 1 column width 85 mm = 3.34 in
## 1.5 column width 114 mm = 4.49 in
## 2 columns width 174 mm = 6.85 in
## max height: 11 in

pdf(file=paste(pathFigures, "Figure2.pdf",sep=""), width=6.85, height=6)

m=matrix(rep(NA, 12*15), nrow=12)

for(i in 1:6){
  m[i,]=c(rep(1,4), rep(2,4), rep(3,7))
}

for(i in 7:12){
  m[i,]=c(rep(4,4), rep(5,4), rep(6,7))
}

layout(m)

par(oma=c(0,0,0,1.5))

#####################################################################################

labels=list()
labels[["human"]]=c("A", "B", "C")
labels[["mouse"]]=c("D", "E", "F")
samples=c("heart", "limb")
names(samples)=c("human", "mouse")
ypos=c(0.75,0.25)
names(ypos)=c("human", "mouse")

for(sp in c("human", "mouse")){

  signif.great1Mb=results[[sp]][["signif.great1Mb"]]
  great1Mb=results[[sp]][["great1Mb"]]

  ## barplot - significant categories, GREAT 1Mb

  top10.great=signif.great1Mb[10:1]
  log10fdr=-log10(great1Mb[top10.great,"FDR"])

  par(mar=c(3.1, 0.5, 2.1, 3.5))
  b=barplot(log10fdr, horiz=T, axes=F, col="steelblue", border="steelblue", space=1, xlim=c(0,10))
  axis(side=1, mgp=c(3,0.5,0), cex.axis=0.9)

  mtext("-log10 FDR", side=1, line=1.5, cex=0.6)

  text(great1Mb[top10.great, "GOName"], y=b+1, x=0, adj=c(0,0.5), cex=0.9, xpd=NA)

  mtext(labels[[sp]][1], side=3, at=-0.09, line=0.5, font=2)

  mtext("GREAT", side=3, at=5.5, cex=0.7, line=0.5)

  #####################################################################################

  ## barplot - significant categories, contacts 1Mb

  signif.contacts1Mb=results[[sp]][["signif.contacts1Mb"]]
  contacts1Mb=results[[sp]][["contacts1Mb"]]

  top10.contacts=signif.contacts1Mb[10:1]
  log10fdr=-log10(contacts1Mb[top10.contacts,"FDR"])

  par(mar=c(3.1, 0.5, 2.1, 3.5))
  b=barplot(log10fdr, horiz=T, axes=F, col="steelblue", border="steelblue", space=1, xlim=c(0,10))
  axis(side=1, mgp=c(3,0.5,0), cex.axis=0.9)

  mtext("-log10 FDR", side=1, line=1.5, cex=0.6)

  text(contacts1Mb[top10.contacts, "GOName"], y=b+1, x=0, adj=c(0,0.5), cex=0.9, xpd=NA)

  mtext(labels[[sp]][2], side=3, at=-0.15, line=0.5, font=2)

  mtext("GOntact", side=3, at=5.5, cex=0.7, line=0.5)

  #####################################################################################

  ## comparison between GREAT and GOntact

  signif.any1Mb=results[[sp]][["signif.any1Mb"]]

  par(mar=c(3.1, 5.0,2.1,0.5))

  col.signifany=rep("black", length(signif.any1Mb))
  names(col.signifany)=signif.any1Mb

  col.signifany[signif.great1Mb]="red"
  col.signifany[signif.contacts1Mb]="steelblue"
  col.signifany[intersect(signif.great1Mb, signif.contacts1Mb)]="darkorange"

  plot(100*great1Mb[signif.any1Mb, "Observed"], 100*contacts1Mb[signif.any1Mb, "Observed"], pch=20, axes=F, xlab="", ylab="", main="", xlim=c(0,25), ylim=c(0,25), col=col.signifany)

  abline(0, 1, col="black", lty=3)

  axis(side=1, mgp=c(3,0.5,0), cex.axis=0.9)
  axis(side=2, mgp=c(3,0.75,0), cex.axis=0.9)
  mtext("% enh. associated w. category, GREAT", side=1, line=1.75, cex=0.7)
  mtext("% enh. associated w. category, GOntact", side=2, line=2.25, cex=0.7)

  box()

  if(sp=="human"){
    legend("topleft", col=c("red", "steelblue", "darkorange"), pch=20, legend=c("GREAT FDR<0.01", "GOntact FDR<0.01", "both methods FDR<0.01"), inset=0.01)
  }

  mtext(labels[[sp]][3], side=3, at=-5.15, line=0.5, font=2)

  ## label for species and sample

  mtext(paste(sp, samples[sp]), side=4, line=0.15, outer=T, at=ypos[sp], cex=0.7)

  #####################################################################################
}

dev.off()

#####################################################################################
