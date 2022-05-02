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
    great=enrichment.results[[sp]][[sample]][[domain]][["GREAT_upstream5kb_downstream1kb_extend1Mb"]]
    rownames(great)=great$GOID
    great$Enrichment=great$Observed/great$Expected

    contacts=enrichment.results[[sp]][[sample]][[domain]][["contacts_mindist0kb_maxdist1Mb"]]
    rownames(contacts)=contacts$GOID
    contacts$Enrichment=contacts$Observed/contacts$Expected

    hybrid=enrichment.results[[sp]][[sample]][[domain]][["hybrid_mindist25kb_maxdist1Mb"]]
    rownames(hybrid)=hybrid$GOID
    hybrid$Enrichment=hybrid$Observed/hybrid$Expected

    ## replace very small FDRs

    great$FDR[which(great$FDR<minFDR)]=minFDR
    contacts$FDR[which(contacts$FDR<minFDR)]=minFDR
    hybrid$FDR[which(hybrid$FDR<minFDR)]=minFDR

    ## signif categories
    signif.great=great$GOID[which(great$FDR<maxFDR  & great$Enrichment>=minEnrichment)]
    signif.contacts=contacts$GOID[which(contacts$FDR<maxFDR & contacts$Enrichment>=minEnrichment)]
    signif.hybrid=hybrid$GOID[which(hybrid$FDR<maxFDR & hybrid$Enrichment>=minEnrichment)]

    ## smaller names

    great$GOName[which(great$GOName=="positive regulation of epithelial cell proliferation")]="positive reg. epithelial cell proliferation"
    great$GOName[which(great$GOName=="regulation of cellular macromolecule biosynthetic process")]="macromolecule biosynthetic process"
    contacts$GOName[which(contacts$GOName=="positive regulation of transcription by RNA polymerase II")]="positive reg. of transcription by RNA pol II"

    ## save results

    results[[sp]]=list("great"=great, "contacts"=contacts, "hybrid"=hybrid, "signif.great"=signif.great, "signif.contacts"=signif.contacts, "signif.hybrid"=signif.hybrid)

  }

  prepare=FALSE
}

#####################################################################################

## 1 column width 85 mm = 3.34 in
## 1.5 column width 114 mm = 4.49 in
## 2 columns width 174 mm = 6.85 in
## max height: 11 in

pdf(file=paste(pathFigures, "Figure2.pdf",sep=""), width=6.85, height=9)

m=matrix(rep(NA, 18*16), nrow=18)

for(i in 1:6){
  m[i,]=c(rep(1,4), rep(2,4), rep(7,4), rep(8,4))
}

for(i in 7:12){
  m[i,]=c(rep(3,4), rep(4,4), rep(9,4), rep(10,4))
}

for(i in 13:18){
  m[i,]=c(rep(5,4), rep(6,4), rep(11,4), rep(12,4))
}

layout(m)

par(oma=c(0,0,0,1.5))

#####################################################################################

labels=list()
labels[["human"]]=c("A", "C", "E")
labels[["mouse"]]=c("B", "D", "F")
samples=c("heart", "heart")
names(samples)=c("human", "mouse")
ypos=c(0.75,0.25)
names(ypos)=c("human", "mouse")

prop.xax=seq(from=-25, to=0, by=5)
prop.xax.lab=-prop.xax
prop.xlim=c(-27,0)

for(sp in c("human", "mouse")){

  ## significant categories, GREAT

  signif.great=results[[sp]][["signif.great"]]
  great=results[[sp]][["great"]]
  top10.great=signif.great[10:1]

  ## barplot, observed vs expected
  observed=great[top10.great,"Observed"]
  expected=great[top10.great,"Expected"]

  all.values=-100*as.numeric(rbind(expected, observed))

  par(mar=c(4.1, 1.5, 2.1, 0.35))
  b=barplot(all.values, horiz=T, axes=F, col=c("black","indianred"), border=NA, space=c(2.25,0.35), xlim=prop.xlim)

  axis(side=1, at=prop.xax, labels=prop.xax.lab, mgp=c(3,0.5,0), cex.axis=0.9)

  mtext(labels[[sp]][1], side=3, at=prop.xlim[1], line=0.5, font=2)
  mtext("% associated enhancers", side=1, line=1.75, cex=0.65)

  if(sp=="human"){
    legend("bottomleft", fill=c("indianred", "black"), border=c("indianred", "black"), legend=c("observed", "expected"), inset=c(-0.175, 0.02),bty="n", cex=0.95, xpd=NA)
  }


  ## barplot, log10 FDR

  log10fdr=-log10(great[top10.great,"FDR"])

  par(mar=c(4.1, 0.35, 2.1, 1.5))
  b=barplot(log10fdr, horiz=T, axes=F, col="steelblue", border=NA, space=1, xlim=c(0,10))

  xax=seq(from=0, to=10, by=2)
  xax.labels=as.character(xax)
  xax.labels[length(xax.labels)]=">10"
  axis(side=1, mgp=c(3,0.5,0), cex.axis=0.9, at=xax, labels=xax.labels)

  mtext("-log10 FDR", side=1, line=1.75, cex=0.65)

  text(great[top10.great, "GOName"], y=b+1, x=0, adj=c(0.5,0.5), cex=0.92, xpd=NA)


  mtext(paste("GREAT,", sp, samples[sp]) , side=3, at=0, cex=0.7, line=0.5)

  #####################################################################################

  ## significant categories, contacts

  signif.contacts=results[[sp]][["signif.contacts"]]
  contacts=results[[sp]][["contacts"]]

  top10.contacts=signif.contacts[10:1]

  ## barplot, observed vs expected
  observed=contacts[top10.contacts,"Observed"]
  expected=contacts[top10.contacts,"Expected"]

  all.values=-100*as.numeric(rbind(expected, observed))

  par(mar=c(4.1, 1.5, 2.1, 0.35))
  b=barplot(all.values, horiz=T, axes=F, col=c("black", "indianred"), border=NA, space=c(2.25, 0.35), xlim=prop.xlim)
  axis(side=1, at=prop.xax, labels=prop.xax.lab, mgp=c(3,0.5,0), cex.axis=0.9)

  mtext(labels[[sp]][2], side=3, at=prop.xlim[1], line=0.5, font=2)
  mtext("% associated enhancers", side=1, line=1.75, cex=0.65)

 ## barplot, log10 FDR

  log10fdr=-log10(contacts[top10.contacts,"FDR"])

  par(mar=c(4.1, 0.35, 2.1, 1.5))
  b=barplot(log10fdr, horiz=T, axes=F, col="steelblue", border="steelblue", space=1, xlim=c(0,10))

  xax=seq(from=0, to=10, by=2)
  xax.labels=as.character(xax)
  xax.labels[length(xax.labels)]=">10"
  axis(side=1, mgp=c(3,0.5,0), cex.axis=0.9, at=xax, labels=xax.labels)

  mtext("-log10 FDR", side=1, line=1.75, cex=0.65)

  text(contacts[top10.contacts, "GOName"], y=b+1, x=0, adj=c(0.5,0.5), cex=0.92, xpd=NA)

  mtext(paste("GOntact,", sp, samples[sp]), side=3, at=0, cex=0.7, line=0.5)

#####################################################################################

  ## significant categories, hybrid

  signif.hybrid=results[[sp]][["signif.hybrid"]]
  hybrid=results[[sp]][["hybrid"]]

  top10.hybrid=signif.hybrid[10:1]

  ## barplot, observed vs expected
  observed=hybrid[top10.hybrid,"Observed"]
  expected=hybrid[top10.hybrid,"Expected"]

  all.values=-100*as.numeric(rbind(expected, observed))


  par(mar=c(4.1, 1.5, 2.1, 0.35))
  b=barplot(all.values, horiz=T, axes=F, col=c("black", "indianred"), border=NA, space=c(2.25, 0.35), xlim=prop.xlim)
  axis(side=1, at=prop.xax, labels=prop.xax.lab, mgp=c(3,0.5,0), cex.axis=0.9)

  mtext(labels[[sp]][3], side=3, at=prop.xlim[1], line=0.5, font=2)
  mtext("% associated enhancers", side=1, line=1.75, cex=0.65)

 ## barplot, log10 FDR

  log10fdr=-log10(hybrid[top10.hybrid,"FDR"])

  par(mar=c(4.1, 0.35, 2.1, 1.5))
  b=barplot(log10fdr, horiz=T, axes=F, col="steelblue", border="steelblue", space=1, xlim=c(0,10))

  xax=seq(from=0, to=10, by=2)
  xax.labels=as.character(xax)
  xax.labels[length(xax.labels)]=">10"
  axis(side=1, mgp=c(3,0.5,0), cex.axis=0.9, at=xax, labels=xax.labels)

  mtext("-log10 FDR", side=1, line=1.75, cex=0.65)

  text(hybrid[top10.hybrid, "GOName"], y=b+1, x=0, adj=c(0.5,0.5), cex=0.92, xpd=NA)

  mtext(paste("hybrid,", sp, samples[sp]), side=3, at=0, cex=0.7, line=0.5)


  #####################################################################################
}

dev.off()

#####################################################################################
