####################################################################################

pathFigures="../../data_for_publication/figures/"

options(stringsAsFactors=F)

#####################################################################################

if(load){
  load(paste(pathFigures, "RData/data.enrichment.results.RData",sep=""))
  load(paste(pathFigures, "RData/data.foreground.GO.association.RData",sep=""))

  species=c("human", "mouse")
  samples=c("Vista_heart_vs_ENCODE", "Vista_limb_vs_ENCODE")
  domain="biological_process"

  minFDR=1e-10
  maxFDR=0.01
  minEnrichment=2

  load=FALSE
}

#####################################################################################

if(prepare){

  results=list()

  for(i in 1:2){
    sp=species[i]
    sample=samples[i]

    test.great1Mb=enrichment.results[[sp]][[sample]][[domain]][["GREAT_upstream5kb_downstream1kb_extend1Mb"]]
    test.great1Mb$Enrichment=test.great1Mb$Observed/test.great1Mb$Expected
    rownames(test.great1Mb)=test.great1Mb$GOID

    assoc.great1Mb=foreground.association[[sp]][[sample]][[domain]][["GREAT_upstream5kb_downstream1kb_extend1Mb"]]

    test.contacts1Mb=enrichment.results[[sp]][[sample]][[domain]][["contacts_mindist0kb_maxdist1Mb"]]
    test.contacts1Mb$Enrichment=test.contacts1Mb$Observed/test.contacts1Mb$Expected
    rownames(test.contacts1Mb)=test.contacts1Mb$GOID

    assoc.contacts1Mb=foreground.association[[sp]][[sample]][[domain]][["contacts_mindist0kb_maxdist1Mb"]]

    ## significant results

    signif.great1Mb=test.great1Mb$GOID[which(test.great1Mb$FDR<maxFDR & test.great1Mb$Enrichment>=minEnrichment)]
    signif.contacts1Mb=test.contacts1Mb$GOID[which(test.contacts1Mb$FDR<maxFDR & test.contacts1Mb$Enrichment>=minEnrichment)]
    signif.both=intersect(signif.great1Mb, signif.contacts1Mb)

    nb.enh.great=test.great1Mb[signif.both, "NbElementsForeground"]
    nb.enh.contacts=test.contacts1Mb[signif.both, "NbElementsForeground"]
    nb.enh.common=unlist(lapply(signif.both, function(x) {y=unlist(strsplit(assoc.contacts1Mb[which(assoc.contacts1Mb$V1==x),2], split=",")) ; z=unlist(strsplit(assoc.great1Mb[which(assoc.great1Mb$V1==x),2], split=",")) ; return(length(intersect(y,z))) }))
    nb.enh.any=unlist(lapply(signif.both, function(x) {y=unlist(strsplit(assoc.contacts1Mb[which(assoc.contacts1Mb$V1==x),2], split=",")) ; z=unlist(strsplit(assoc.great1Mb[which(assoc.great1Mb$V1==x),2], split=",")) ; return(length(unique(c(y,z))))}))
   }

  prepare=FALSE
}

#####################################################################################



#####################################################################################
