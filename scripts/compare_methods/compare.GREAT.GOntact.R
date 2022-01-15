##################################################################################

pathGREAT="../../results/GO_enrichment_GREAT/"
pathGOntact="../../results/GO_enrichment_contacts/"
pathFigures="../../results/figures/"

##################################################################################

sp="human"
dataset="FANTOM5.first1000.kidney.enhancers.hg38"
space="biological_process"

maxFDR=0.01

##################################################################################

great=read.table(paste(pathGREAT, sp,"/classical_upstream5kb_downstream1kb_extend1Mb/",dataset,"/enrichment_",space,"_Ensembl94_backgroundFANTOM5.Laverre2022.txt",sep=""), h=T, stringsAsFactors=F, sep="\t", quote="\"")

gontact=read.table(paste(pathGOntact, sp, "/mindistOKb_maxdist1Mb_extendOverlap5Kb/FANTOM5.first1000.kidney.enhancers.hg38/GOEnrichment/uniqueGO.biological_process.Background.EnhancersContactsFANTOM5.Laverre2022.EnhancerScale.EnhancerCountOnce.txt", sep=""), h=T, stringsAsFactors=F, sep="\t", quote="\"")

##################################################################################

length(which(great$FDR<maxFDR))
length(which(gontact$FDR<maxFDR))

##################################################################################

great2=read.table(paste(pathGREAT, sp,"/classical_mindist25kb_extend1Mb/",dataset,"/enrichment_",space,"_Ensembl94_backgroundFANTOM5.Laverre2022.txt",sep=""), h=T, stringsAsFactors=F, sep="\t", quote="\"")

gontact2=read.table(paste(pathGOntact, sp, "/mindist25Kb_maxdist1Mb_extendOverlap5Kb/FANTOM5.first1000.kidney.enhancers.hg38/GOEnrichment/uniqueGO.biological_process.Background.EnhancersContactsFANTOM5.Laverre2022.EnhancerScale.EnhancerCountOnce.txt", sep=""), h=T, stringsAsFactors=F, sep="\t", quote="\"")

length(which(great2$FDR<maxFDR))
length(which(gontact2$FDR<maxFDR))

##################################################################################
