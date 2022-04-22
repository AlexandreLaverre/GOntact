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

  prepare=TRUE
}

#####################################################################################



#####################################################################################
