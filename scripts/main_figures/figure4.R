#####################################################################################

pathFigures="../../data_for_publication/figures/"
pathJASPAR="../../data/JASPAR/"

options(stringsAsFactors=F)

library(imager)

#####################################################################################

if(load){

  load(paste(pathFigures, "RData/data.motif.enrichment.RData",sep=""))

  load=FALSE
}

#####################################################################################
