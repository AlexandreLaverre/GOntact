#####################################################################################

pathResults="../../results/"
pathFigures="../../data_for_publication/figures/"

options(stringsAsFactors=F)

#####################################################################################

enrichment.results=list()

for(sp in c("human", "mouse")){
  enrichment.results[[sp]]=list()

  datasets=system(paste("ls ", pathResults, sp,sep=""), intern=T)

  for(data in datasets){
    enrichment.results[[sp]][[data]]=list()

    domains=system(paste("ls ", pathResults, sp, "/", data, sep=""), intern=T)

    for(domain in domains){
      enrichment.results[[sp]][[data]][[domain]]=list()

      methods=system(paste("ls ", pathResults, sp, "/", data, "/", domain, " | grep -v txt", sep=""), intern=T)

      for(method in methods){
        res=read.table(paste(pathResults, sp, "/", data, "/", domain, "/", method, "/GOntact_enrichment_results.txt",sep=""), h=T,sep="\t", quote="\"")

        enrichment.results[[sp]][[data]][[domain]][[method]]=res
      }
    }
  }
}

#####################################################################################

save(enrichment.results, file=paste(pathFigures, "RData/data.enrichment.results.RData",sep=""))

#####################################################################################
