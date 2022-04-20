#####################################################################################

pathFigures="../../data_for_publication/figures/"

options(stringsAsFactors=F)

#####################################################################################

if(load){
  load(paste(pathFigures, "RData/data.enrichment.results.RData",sep=""))

  sp="human"
  sample="Vista_heart_vs_ENCODE"
  domain="biological_process"

  maxFDR=0.01
  minEnrichment=2

  load=FALSE
}

#####################################################################################

if(prepare){

  ## enrichment results
  great1Mb=enrichment.results[[sp]][[sample]][[domain]][["GREAT_upstream5kb_downstream1kb_extend1Mb"]]
  rownames(great1Mb)=great1Mb$GOID
  great1Mb$Enrichment=great1Mb$Observed/great1Mb$Expected

  great2Mb=enrichment.results[[sp]][[sample]][[domain]][["GREAT_upstream5kb_downstream1kb_extend2Mb"]]
  rownames(great2Mb)=great2Mb$GOID
  great2Mb$Enrichment=great2Mb$Observed/great2Mb$Expected

  contacts1Mb=enrichment.results[[sp]][[sample]][[domain]][["contacts_mindist25kb_maxdist1Mb"]]
  rownames(contacts1Mb)=contacts1Mb$GOID
  contacts1Mb$Enrichment=contacts1Mb$Observed/contacts1Mb$Expected

  contacts2Mb=enrichment.results[[sp]][[sample]][[domain]][["contacts_mindist25kb_maxdist2Mb"]]
  rownames(contacts2Mb)=contacts2Mb$GOID
  contacts2Mb$Enrichment=contacts2Mb$Observed/contacts2Mb$Expected

  hybrid1Mb=enrichment.results[[sp]][[sample]][[domain]][["hybrid_mindist25kb_maxdist1Mb"]]
  rownames(hybrid1Mb)=hybrid1Mb$GOID
  hybrid1Mb$Enrichment=hybrid1Mb$Observed/hybrid1Mb$Expected

  ## signif categories
  signif.great1Mb=great1Mb$GOID[which(great1Mb$FDR<maxFDR & great1Mb$Enrichment>=minEnrichment)]
  signif.great2Mb=great2Mb$GOID[which(great2Mb$FDR<maxFDR & great2Mb$Enrichment>=minEnrichment)]
  signif.contacts1Mb=contacts1Mb$GOID[which(contacts1Mb$FDR<maxFDR & contacts1Mb$Enrichment>=minEnrichment)]
  signif.contacts2Mb=contacts2Mb$GOID[which(contacts2Mb$FDR<maxFDR & contacts2Mb$Enrichment>=minEnrichment)]
  signif.hybrid1Mb=hybrid1Mb$GOID[which(hybrid1Mb$FDR<maxFDR & hybrid1Mb$Enrichment>=minEnrichment)]


  ## direct comparison great vs contacts, 1Mb
  common1Mb=intersect(great1Mb$GOID, contacts1Mb$GOID)

  load=FALSE
}

#####################################################################################

## 1 column width 85 mm = 3.34 in
## 1.5 column width 114 mm = 4.49 in
## 2 columns width 174 mm = 6.85 in
## max height: 11 in

pdf(file=paste(pathFigures, "Figure1.pdf",sep=""), width=6.85, height=5.5)

m=matrix(rep(NA, 12*9), nrow=12)

for(i in 1:5){
  m[i,]=rep(1,9)
}

for(i in 6:12){
  m[i,]=c(rep(2,2), rep(3,2), rep(4,5))
}

layout(m)

#####################################################################################



#####################################################################################



#####################################################################################

dev.off()

#####################################################################################
