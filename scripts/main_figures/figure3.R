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

  diff.great.contacts=list()
  prop.common=list()

  id=c()

  nb.signif.great=c()
  nb.signif.contacts=c()
  nb.signif.common=c()
  nb.signif.any=c()

  for(sp in species){
    for(sample in names(enrichment.results[[sp]])){
      test.great1Mb=enrichment.results[[sp]][[sample]][[domain]][["GREAT_upstream5kb_downstream1kb_extend1Mb"]]
      test.great1Mb$Enrichment=test.great1Mb$Observed/test.great1Mb$Expected
      rownames(test.great1Mb)=test.great1Mb$GOID

      test.contacts1Mb=enrichment.results[[sp]][[sample]][[domain]][["contacts_mindist0kb_maxdist1Mb"]]
      test.contacts1Mb$Enrichment=test.contacts1Mb$Observed/test.contacts1Mb$Expected
      rownames(test.contacts1Mb)=test.contacts1Mb$GOID

      ## signif categories
      signif.great1Mb=test.great1Mb$GOID[which(test.great1Mb$FDR<maxFDR & test.great1Mb$Enrichment>=minEnrichment)]
      signif.contacts1Mb=test.contacts1Mb$GOID[which(test.contacts1Mb$FDR<maxFDR & test.contacts1Mb$Enrichment>=minEnrichment)]
      signif.any1Mb=unique(c(signif.great1Mb, signif.contacts1Mb))
      signif.both1Mb=intersect(signif.great1Mb, signif.contacts1Mb)

      nb.signif.great=c(nb.signif.great, length(signif.great1Mb))
      nb.signif.contacts=c(nb.signif.contacts, length(signif.contacts1Mb))
      nb.signif.common=c(nb.signif.common, length(signif.both1Mb))
      nb.signif.any=c(nb.signif.any, length(signif.any1Mb))

      id=c(id, paste(sp, sample))


      ## assoc.great1Mb=foreground.association[[sp]][[sample]][[domain]][["GREAT_upstream5kb_downstream1kb_extend1Mb"]]
      ## assoc.contacts1Mb=foreground.association[[sp]][[sample]][[domain]][["contacts_mindist0kb_maxdist1Mb"]]

      ## ## significant results

      ## nb.enh.contacts=test.contacts1Mb[signif.both, "NbElementsForeground"]
      ## nb.enh.common=unlist(lapply(signif.both, function(x) {y=unlist(strsplit(assoc.contacts1Mb[which(assoc.contacts1Mb$V1==x),2], split=",")) ; z=unlist(strsplit(assoc.great1Mb[which(assoc.great1Mb$V1==x),2], split=",")) ; return(length(intersect(y,z))) }))
      ## nb.enh.any=unlist(lapply(signif.both, function(x) {y=unlist(strsplit(assoc.contacts1Mb[which(assoc.contacts1Mb$V1==x),2], split=",")) ; z=unlist(strsplit(assoc.great1Mb[which(assoc.great1Mb$V1==x),2], split=",")) ; return(length(unique(c(y,z))))}))

      ## diff.great.contacts[[sp]][[sample]]=nb.enh.great-nb.enh.contacts
      ## prop.common[[sp]][[sample]]=nb.enh.common/nb.enh.any
    }
  }

  prepare=FALSE
}

#####################################################################################

## 1 column width 85 mm = 3.34 in
## 1.5 column width 114 mm = 4.49 in
## 2 columns width 174 mm = 6.85 in
## max height: 11 in

pdf(file=paste(pathFigures, "Figure3.pdf",sep=""), width=6.85, height=6)

#####################################################################################

## ## boxplot diff GREAT - contacts

## nbsamples=length(diff.great.contacts[["human"]])+length(diff.great.contacts[["mouse"]])
## xlim=c(0.5, nbsamples+0.5)

## all.values=c(unlist(diff.great.contacts[["human"]]), unlist(diff.great.contacts[["mouse"]]))
## ylim=range(all.values)

## plot(1, xlim=xlim, ylim=ylim, type="n", axes=F, xlab="", ylab="", main="")

## xpos=1

## for(sample in names(diff.great.contacts[["human"]])){
##   this.val=diff.great.contacts[["human"]][[sample]]
##   boxplot(this.val, at=xpos, add=T, axes=F)
##   xpos=xpos+1
## }

## for(sample in names(diff.great.contacts[["mouse"]])){
##   this.val=diff.great.contacts[["mouse"]][[sample]]
##   boxplot(this.val, at=xpos, add=T, axes=F)
##   xpos=xpos+1
## }

## axis(side=1, at=1:nbsamples, labels=rep("", nbsamples))
## axis(side=2)

#####################################################################################

dev.off()

#####################################################################################
