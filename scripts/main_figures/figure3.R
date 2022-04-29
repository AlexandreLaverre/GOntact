####################################################################################

pathFigures="../../data_for_publication/figures/"

options(stringsAsFactors=F)

#####################################################################################

if(load){
  load(paste(pathFigures, "RData/data.enrichment.results.RData",sep=""))
  load(paste(pathFigures, "RData/data.foreground.GO.association.RData",sep=""))

  species=c("human", "mouse")
  domain="biological_process"

  minFDR=1e-10
  maxFDR=0.05
  minEnrichment=1.5

  load=FALSE
}

#####################################################################################

if(prepare){

  id=c()

  nb.signif.great1Mb=c()
  nb.signif.contacts1Mb=c()
  nb.signif.common1Mb=c()
  nb.signif.any1Mb=c()

  nb.signif.great2Mb=c()
  nb.signif.contacts2Mb=c()
  nb.signif.common2Mb=c()
  nb.signif.any2Mb=c()

  nb.great1Mb.common=list()
  nb.contacts1Mb.common=list()
  nb.shared1Mb.common=list()

  for(sp in species){
    for(sample in names(enrichment.results[[sp]])){
      ## 1Mb
      test.great1Mb=enrichment.results[[sp]][[sample]][[domain]][["GREAT_upstream5kb_downstream1kb_extend1Mb"]]
      test.great1Mb$Enrichment=test.great1Mb$Observed/test.great1Mb$Expected
      rownames(test.great1Mb)=test.great1Mb$GOID

      test.contacts1Mb=enrichment.results[[sp]][[sample]][[domain]][["contacts_mindist0kb_maxdist1Mb"]]
      test.contacts1Mb$Enrichment=test.contacts1Mb$Observed/test.contacts1Mb$Expected
      rownames(test.contacts1Mb)=test.contacts1Mb$GOID

      ## 2Mb
      test.great2Mb=enrichment.results[[sp]][[sample]][[domain]][["GREAT_upstream5kb_downstream1kb_extend2Mb"]]
      test.great2Mb$Enrichment=test.great2Mb$Observed/test.great2Mb$Expected
      rownames(test.great2Mb)=test.great2Mb$GOID

      test.contacts2Mb=enrichment.results[[sp]][[sample]][[domain]][["contacts_mindist0kb_maxdist2Mb"]]
      test.contacts2Mb$Enrichment=test.contacts2Mb$Observed/test.contacts2Mb$Expected
      rownames(test.contacts2Mb)=test.contacts2Mb$GOID

      ## signif categories

      ## 1Mb
      signif.great1Mb=test.great1Mb$GOID[which(test.great1Mb$FDR<maxFDR & test.great1Mb$Enrichment>=minEnrichment)]
      signif.contacts1Mb=test.contacts1Mb$GOID[which(test.contacts1Mb$FDR<maxFDR & test.contacts1Mb$Enrichment>=minEnrichment)]
      signif.any1Mb=unique(c(signif.great1Mb, signif.contacts1Mb))
      signif.both1Mb=intersect(signif.great1Mb, signif.contacts1Mb)

      nb.signif.great1Mb=c(nb.signif.great1Mb, length(signif.great1Mb))
      nb.signif.contacts1Mb=c(nb.signif.contacts1Mb, length(signif.contacts1Mb))
      nb.signif.common1Mb=c(nb.signif.common1Mb, length(signif.both1Mb))
      nb.signif.any1Mb=c(nb.signif.any1Mb, length(signif.any1Mb))

      ## 2Mb
      signif.great2Mb=test.great2Mb$GOID[which(test.great2Mb$FDR<maxFDR & test.great2Mb$Enrichment>=minEnrichment)]
      signif.contacts2Mb=test.contacts2Mb$GOID[which(test.contacts2Mb$FDR<maxFDR & test.contacts2Mb$Enrichment>=minEnrichment)]
      signif.any2Mb=unique(c(signif.great2Mb, signif.contacts2Mb))
      signif.both2Mb=intersect(signif.great2Mb, signif.contacts2Mb)

      nb.signif.great2Mb=c(nb.signif.great2Mb, length(signif.great2Mb))
      nb.signif.contacts2Mb=c(nb.signif.contacts2Mb, length(signif.contacts2Mb))
      nb.signif.common2Mb=c(nb.signif.common2Mb, length(signif.both2Mb))
      nb.signif.any2Mb=c(nb.signif.any2Mb, length(signif.any2Mb))

      id=c(id, paste(sp, sample))

      assoc.great1Mb=foreground.association[[sp]][[sample]][[domain]][["GREAT_upstream5kb_downstream1kb_extend1Mb"]]
      assoc.contacts1Mb=foreground.association[[sp]][[sample]][[domain]][["contacts_mindist0kb_maxdist1Mb"]]

      nb.enh.great=test.great1Mb[signif.both1Mb, "NbElementsForeground"]
      nb.enh.contacts=test.contacts1Mb[signif.both1Mb, "NbElementsForeground"]
      nb.enh.common=unlist(lapply(signif.both1Mb, function(x) {y=unlist(strsplit(assoc.contacts1Mb[which(assoc.contacts1Mb$V1==x),2], split=",")) ; z=unlist(strsplit(assoc.great1Mb[which(assoc.great1Mb$V1==x),2], split=",")) ; return(length(intersect(y,z))) }))
      nb.enh.any=nb.enh.great+nb.enh.contacts-nb.enh.common

      names(nb.enh.great)=signif.both1Mb
      names(nb.enh.contacts)=signif.both1Mb
      names(nb.enh.common)=signif.both1Mb
      names(nb.enh.any)=signif.both1Mb

      nb.great1Mb.common[[paste(sp, sample)]]=nb.enh.great
      nb.contacts1Mb.common[[paste(sp, sample)]]=nb.enh.contacts
      nb.shared1Mb.common[[paste(sp, sample)]]=nb.enh.common
    }
  }

  prepare=FALSE
}

#####################################################################################

## 1 column width 85 mm = 3.34 in
## 1.5 column width 114 mm = 4.49 in
## 2 columns width 174 mm = 6.85 in
## max height: 11 in

pdf(file=paste(pathFigures, "Figure3.pdf",sep=""), width=4.95, height=6)

m=matrix(rep(NA, 2*9), nrow=2)

m[1,]=c(rep(1,3), rep(2, 3), rep(3,3))
m[2,]=c(rep(4,3), rep(5, 3), rep(6,3))

layout(m)

#####################################################################################

sample="Vista_heart_vs_ENCODE"

labels=list("human"=c("A", "B", "C"), "mouse"=c("D", "E", "F"))

for(sp in c("human", "mouse")){
  this.id=paste(sp, sample)

  ## barplot nb great, contacts, shared

  nbg=nb.signif.great1Mb[which(id==this.id)]
  nbc=nb.signif.contacts1Mb[which(id==this.id)]
  nbs=nb.signif.common1Mb[which(id==this.id)]

  par(mar=c(5.1, 4.1, 2.5, 1.5))
  b=barplot(c(nbg, nbc, nbs), col=c("black", "darkorange", "navy"), density=20, border=c("black", "darkorange", "navy"), xlab="", names=rep("", 3), space=c(1, 0.5, 0.5))

  mtext("nb. significant GO categories", side=2, line=2.5, cex=0.8)

  mtext(c("GREAT", "GOntact", "shared"),side=1, at=b, las=2,line=0.1, cex=0.75)

  mtext(labels[[sp]][1], side=3, at=-0.95, font=2, line=1, cex=1.1)


  ## boxplot nb contacted enhancers


  par(mar=c(5.1, 3.5, 2.5, 2.0))

  boxplot(nb.great1Mb.common[[this.id]], nb.contacts1Mb.common[[this.id]], col="white", border=c("black", "darkorange"), boxwex=0.6, names=rep("",2), notch=T, pch=20, axes=F)
  axis(side=1, at=c(1,2), labels=rep("",2))
  axis(side=2)

  mtext("nb. associated enhancers", side=2, line=2.5, cex=0.8)
  mtext(c("GREAT", "GOntact"), side=1, at=c(1,2), line=1, cex=0.75, las=2)

  mtext(labels[[sp]][2], side=3, at=-0.5, font=2, line=1, cex=1.1)

  ## boxplot proportion shared

  par(mar=c(5.1, 3.5, 2.5, 2))

  prop.shared.gontact=nb.shared1Mb.common[[this.id]]/nb.contacts1Mb.common[[this.id]]
  prop.shared.great=nb.shared1Mb.common[[this.id]]/nb.great1Mb.common[[this.id]]

  boxplot(prop.shared.great, prop.shared.gontact, border=c("black", "darkorange"), col="white", notch=T, axes=F, boxwex=0.6)
  axis(side=1, at=c(1, 2), labels=rep("",2))
  axis(side=2)
  mtext("% shared associated enhancers", side=2, line=2.5, cex=0.8)
  mtext(c("GREAT", "GOntact"), side=1, at=c(1,2), line=1, cex=0.75, las=2)

  mtext(labels[[sp]][3], side=3, at=-0.5, font=2, line=1, cex=1.1)

  mtext(paste(sp, "heart"), side=4, line=0.5, cex=0.75)
}

#####################################################################################

dev.off()

#####################################################################################
