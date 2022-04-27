#####################################################################################

pathFigures="../../data_for_publication/figures/"

options(stringsAsFactors=F)

#####################################################################################

if(load){
  load(paste(pathFigures, "RData/data.enrichment.results.RData",sep=""))

  species=c("human", "mouse")
  samples=c("Vista_heart_vs_ENCODE", "Vista_heart_vs_ENCODE")
  domain="biological_process"

  minFDR=1e-10
  maxFDR=0.05
  minEnrichment=1.5

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

    ## replace very small FDRs

    great1Mb$FDR[which(great1Mb$FDR<minFDR)]=minFDR
    contacts1Mb$FDR[which(contacts1Mb$FDR<minFDR)]=minFDR

    ## signif categories
    signif.great1Mb=great1Mb$GOID[which(great1Mb$FDR<maxFDR  & great1Mb$Enrichment>=minEnrichment)]
    signif.contacts1Mb=contacts1Mb$GOID[which(contacts1Mb$FDR<maxFDR & contacts1Mb$Enrichment>=minEnrichment)]
    signif.any1Mb=unique(c(signif.great1Mb, signif.contacts1Mb))

    ## direct comparison great vs contacts, 1Mb
    common1Mb=intersect(great1Mb$GOID, contacts1Mb$GOID)

    ## smaller names

    great1Mb$GOName[which(great1Mb$GOName=="positive regulation of epithelial cell proliferation")]="positive reg. epithelial cell proliferation"
    great1Mb$GOName[which(great1Mb$GOName=="regulation of cellular macromolecule biosynthetic process")]="macromolecule biosynthetic process"
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

pdf(file=paste(pathFigures, "Figure2.pdf",sep=""), width=6.85, height=7)

m=matrix(rep(NA, 12*16), nrow=12)

for(i in 1:6){
  m[i,]=c(rep(1,4), rep(2,4), rep(3,4), rep(4,4))
}

for(i in 7:12){
  m[i,]=c(rep(5,4), rep(6,4), rep(7,4), rep(8,4))
}

layout(m)

par(oma=c(0,0,0,1.5))

#####################################################################################

labels=list()
labels[["human"]]=c("A", "B")
labels[["mouse"]]=c("C", "D")
samples=c("heart", "heart")
names(samples)=c("human", "mouse")
ypos=c(0.75,0.25)
names(ypos)=c("human", "mouse")

for(sp in c("human", "mouse")){

  ## significant categories, GREAT 1Mb

  signif.great1Mb=results[[sp]][["signif.great1Mb"]]
  great1Mb=results[[sp]][["great1Mb"]]
  top10.great=signif.great1Mb[10:1]

  ## barplot, observed vs expected
  observed=great1Mb[top10.great,"Observed"]
  expected=great1Mb[top10.great,"Expected"]

  all.values=-100*as.numeric(rbind(expected, observed))
  xax=pretty(all.values)
  xax.lab=-xax
  xlim=range(all.values)
  xlim[1]=min(c(xlim, xax))-5
  xlim[2]=0

  par(mar=c(4.1, 1.5, 2.1, 0.35))
  b=barplot(all.values, horiz=T, axes=F, col=c("black","indianred"), border=NA, space=c(2.25,0.35), xlim=xlim)

  axis(side=1, at=xax, labels=xax.lab, mgp=c(3,0.5,0), cex.axis=0.9)

  mtext(labels[[sp]][1], side=3, at=xlim[1], line=0.5, font=2)
  mtext("% associated enhancers", side=1, line=1.75, cex=0.65)

  if(sp=="human"){
    legend("bottomleft", fill=c("indianred", "black"), border=c("indianred", "black"), legend=c("observed", "expected"), inset=c(-0.175, 0.02),bty="n", cex=0.95, xpd=NA)
  }


  ## barplot, log10 FDR

  log10fdr=-log10(great1Mb[top10.great,"FDR"])

  par(mar=c(4.1, 0.35, 2.1, 1.5))
  b=barplot(log10fdr, horiz=T, axes=F, col="steelblue", border=NA, space=1, xlim=c(0,10))

  xax=seq(from=0, to=10, by=2)
  xax.labels=as.character(xax)
  xax.labels[length(xax.labels)]=">10"
  axis(side=1, mgp=c(3,0.5,0), cex.axis=0.9, at=xax, labels=xax.labels)

  mtext("-log10 FDR", side=1, line=1.75, cex=0.65)

  text(great1Mb[top10.great, "GOName"], y=b+1, x=0, adj=c(0.5,0.5), cex=0.92, xpd=NA)


  mtext(paste("GREAT,", sp, samples[sp]) , side=3, at=0, cex=0.7, line=0.5)

  #####################################################################################

  ## significant categories, contacts 1Mb

  signif.contacts1Mb=results[[sp]][["signif.contacts1Mb"]]
  contacts1Mb=results[[sp]][["contacts1Mb"]]

  top10.contacts=signif.contacts1Mb[10:1]

  ## barplot, observed vs expected
  observed=contacts1Mb[top10.contacts,"Observed"]
  expected=contacts1Mb[top10.contacts,"Expected"]

  all.values=-100*as.numeric(rbind(expected, observed))
  xax=pretty(all.values)
  xax.lab=-xax
  xlim=range(all.values)
  xlim[1]=min(c(xlim, xax))
  xlim[2]=0

  par(mar=c(4.1, 1.5, 2.1, 0.35))
  b=barplot(all.values, horiz=T, axes=F, col=c("black", "indianred"), border=NA, space=c(2.25, 0.35), xlim=xlim)
  axis(side=1, at=xax, labels=xax.lab, mgp=c(3,0.5,0), cex.axis=0.9)

  mtext(labels[[sp]][2], side=3, at=xlim[1], line=0.5, font=2)
  mtext("% associated enhancers", side=1, line=1.75, cex=0.65)

 ## barplot, log10 FDR

  log10fdr=-log10(contacts1Mb[top10.contacts,"FDR"])

  par(mar=c(4.1, 0.35, 2.1, 1.5))
  b=barplot(log10fdr, horiz=T, axes=F, col="steelblue", border="steelblue", space=1, xlim=c(0,10))

  xax=seq(from=0, to=10, by=2)
  xax.labels=as.character(xax)
  xax.labels[length(xax.labels)]=">10"
  axis(side=1, mgp=c(3,0.5,0), cex.axis=0.9, at=xax, labels=xax.labels)

  mtext("-log10 FDR", side=1, line=1.75, cex=0.65)

  text(contacts1Mb[top10.contacts, "GOName"], y=b+1, x=0, adj=c(0.5,0.5), cex=0.92, xpd=NA)

  mtext(paste("GOntact,", sp, samples[sp]), side=3, at=0, cex=0.7, line=0.5)

  #####################################################################################
}

dev.off()

#####################################################################################
