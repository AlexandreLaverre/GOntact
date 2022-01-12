##################################################################################

pathResults="../../results/GO_enrichment_GREAT/"

##################################################################################

sp="human"
dataset="FANTOM5.first1000.kidney.enhancers.hg38"
space="biological_process"
background="regulatory_regions"

method1="classical_upstream5kb_downstream1kb_extend1Mb"
method2="classical_mindist25kb_extend1Mb"

##################################################################################

## GO enrichment tests

res1=read.table(paste(pathResults, sp,"/",method1,"/",dataset,"/enrichment_regulatory_regions_",space,"_Ensembl94.txt",sep=""), h=T, stringsAsFactors=F, sep="\t", quote=F)

res2=read.table(paste(pathResults, sp,"/",method2,"/",dataset,"/enrichment_regulatory_regions_",space,"_Ensembl94.txt",sep=""), h=T, stringsAsFactors=F, sep="\t", quote=F)


##################################################################################
