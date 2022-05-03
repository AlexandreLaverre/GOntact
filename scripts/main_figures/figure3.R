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

  intersections=list()
  enhancers.great=list()
  enhancers.gontact=list()
  enhancers.shared=list()

  for(sp in species){
    for(sample in names(enrichment.results[[sp]])){
      test.great=enrichment.results[[sp]][[sample]][[domain]][["GREAT_upstream5kb_downstream1kb_extend1Mb"]]
      test.great$Enrichment=test.great$Observed/test.great$Expected
      rownames(test.great)=test.great$GOID

      test.contacts=enrichment.results[[sp]][[sample]][[domain]][["contacts_mindist0kb_maxdist1Mb"]]
      test.contacts$Enrichment=test.contacts$Observed/test.contacts$Expected
      rownames(test.contacts)=test.contacts$GOID

      test.hybrid=enrichment.results[[sp]][[sample]][[domain]][["hybrid_mindist25kb_maxdist1Mb"]]
      test.hybrid$Enrichment=test.hybrid$Observed/test.hybrid$Expected
      rownames(test.hybrid)=test.hybrid$GOID

      signif.great=test.great$GOID[which(test.great$FDR < maxFDR & test.great$Enrichment >= minEnrichment)]
      signif.contacts=test.contacts$GOID[which(test.contacts$FDR < maxFDR & test.contacts$Enrichment >= minEnrichment)]
      signif.hybrid=test.hybrid$GOID[which(test.hybrid$FDR < maxFDR & test.hybrid$Enrichment >= minEnrichment)]

      only.great=setdiff(signif.great, c(signif.contacts, signif.hybrid))
      only.contacts=setdiff(signif.contacts, c(signif.great, signif.hybrid))
      only.hybrid=setdiff(signif.hybrid, c(signif.great, signif.contacts))

      great.contacts=setdiff(intersect(signif.contacts, signif.great), signif.hybrid)
      great.hybrid=setdiff(intersect(signif.hybrid, signif.great), signif.contacts)
      contacts.hybrid=setdiff(intersect(signif.hybrid, signif.contacts), signif.great)

      triple=intersect(intersect(signif.contacts, signif.great), signif.hybrid)

      intersections[[paste(sp, sample)]]=list("only.great"=only.great, "only.contacts"=only.contacts, "only.hybrid"=only.hybrid, "great.contacts"=great.contacts, "great.hybrid"=great.hybrid, "contacts.hybrid"=contacts.hybrid, "triple"=triple)

      ## shared GO categories - associations

      shared=intersect(signif.contacts, signif.great)

      ## GO-enhancer association

      assoc.great=foreground.association[[sp]][[sample]][[domain]][["GREAT_upstream5kb_downstream1kb_extend1Mb"]]
      assoc.gontact=foreground.association[[sp]][[sample]][[domain]][["contacts_mindist0kb_maxdist1Mb"]]

      enh.great=lapply(shared, function(x) unlist(strsplit(assoc.great$V2[which(assoc.great$V1==x)], split=",")))
      names(enh.great)=shared

      enh.gontact=lapply(shared, function(x) unlist(strsplit(assoc.gontact$V2[which(assoc.gontact$V1==x)], split=",")))
      names(enh.gontact)=shared

      enh.shared=lapply(shared, function(x) intersect(enh.great[[x]], enh.gontact[[x]]))
      names(enh.shared)=shared

      enhancers.great[[paste(sp, sample)]]=enh.great
      enhancers.gontact[[paste(sp, sample)]]=enh.gontact
      enhancers.shared[[paste(sp, sample)]]=enh.shared

    }
  }

  prepare=FALSE
}

#####################################################################################

## 1 column width 85 mm = 3.34 in
## 1.5 column width 114 mm = 4.49 in
## 2 columns width 174 mm = 6.85 in
## max height: 11 in

pdf(file=paste(pathFigures, "Figure3.pdf",sep=""), width=6.85, height=7.5)

m=matrix(rep(NA, 6*10), nrow=6)

for(i in 1:2){
  m[i,]=c(rep(1,6), rep(3, 4))
}

m[3,]=c(rep(2,6), rep(3, 4))

for(i in 4:5){
  m[i,]=c(rep(4,6), rep(6, 4))
}

m[6,]=c(rep(5,6), rep(6, 4))

layout(m)

#####################################################################################

sample="Vista_heart_vs_ENCODE"

labels=list("human"=c("A", "B"), "mouse"=c("C", "D"))


for(sp in c("human", "mouse")){
  this.sample=paste(sp, sample)
  this.int=intersections[[this.sample]]

  ## barplot for the numbers of GO categories in each intersection class

  numbers=unlist(lapply(this.int, length))
  names(numbers)=names(this.int)

  xpos=1:length(numbers)
  names(xpos)=c("only.great", "only.contacts", "only.hybrid", "great.contacts", "great.hybrid", "contacts.hybrid", "triple")

  ylim=c(0, max(numbers))
  xlim=c(0.25, length(numbers)+0.75)

  par(mar=c(0.5,4.5, 2.5,0.5))
  plot(1, type="n", xlab="", ylab="", xlim=xlim, ylim=ylim, yaxs="i", xaxs="i", axes=F)

  rect(xpos-0.25, 0, xpos+0.25, numbers, col="gray40", border="black")

  axis(side=2, mgp=c(3, 0.65, 0), cex.axis=1.15)

  mtext("nb. significant GO categories", side=2, line=2.75, cex=0.9)

  mtext(labels[[sp]][1], side=3, at=-0.725, font=2, line=1, cex=1.25)

  ## labels for the intersections

  par(mar=c(1,4.5,0.5,0.5))
  plot(1, type="n", xlab="", ylab="", xlim=xlim, ylim=c(0.5,3.5), xaxs="i", yaxs="i", axes=F)
  smallx=diff(xlim)/40
  segments(xpos[1]-smallx, 1:3, xpos[7]+smallx, 1:3, lty=3)
  points(xpos[c(1, 4, 5,7)], rep(3,4), pch=20, cex=2) ## GREAT
  points(xpos[c(2, 4, 6,7)], rep(2,4), pch=20, cex=2) ## GOntact
  points(xpos[c(3, 5, 6,7)], rep(1,4), pch=20, cex=2) ## hybrid

  segments(xpos[4], 2, xpos[4], 3)
  segments(xpos[5], 1, xpos[5], 3)
  segments(xpos[6], 1, xpos[6], 2)
  segments(xpos[7], 1, xpos[7], 3)


  mtext("GREAT", at=3, side=2, las=2, cex=0.9)
  mtext("GOntact", at=2, side=2, las=2, cex=0.9)
  mtext("hybrid", at=1, side=2, las=2, cex=0.9)

  ## number of enhancers found with just one of the approaches

  enh.great=enhancers.great[[paste(sp, sample)]]
  enh.gontact=enhancers.gontact[[paste(sp, sample)]]
  enh.shared=enhancers.shared[[paste(sp, sample)]]

  nb.only.great=unlist(lapply(names(enh.great), function(x) length(setdiff(enh.great[[x]], enh.gontact[[x]]))))
  nb.only.gontact=unlist(lapply(names(enh.gontact), function(x) length(setdiff(enh.gontact[[x]], enh.great[[x]]))))
  nb.shared=unlist(lapply(names(enh.shared), function(x) length(intersect(enh.gontact[[x]], enh.great[[x]]))))

  par(mar=c(5.1,4.5,2.5,1.5))

  boxplot(nb.only.great, nb.only.gontact, nb.shared, boxwex=0.65, col="white", border=c("black", "darkorange", "slateblue"), pch=20, cex=1.1, axes=F)
  axis(side=2, mgp=c(3, 0.65, 0), cex.axis=1.15)
  axis(side=1, at=1:3, labels=rep("",3))

  mtext("nb. associated enhancers", side=2, line=2.75, cex=0.9)

  mtext(c("GREAT", "GOntact"), at=c(1,2), side=1, line=1.1, cex=0.9)
  mtext(c("only", "only"), at=c(1,2), side=1, line=2.5, cex=0.9)
  mtext("shared", at=3, side=1, line=1.8, cex=0.9)


  mtext(labels[[sp]][2], side=3, at=-0.35, font=2, line=1, cex=1.25)


  mtext(paste(sp, "heart"), side=4, line=0.25)

}

#####################################################################################

dev.off()

#####################################################################################
#####################################################################################
