#####################################################################################

pathResults="../../results/"
pathFigures="../../data_for_publication/figures/"

options(stringsAsFactors=F)

#####################################################################################

regulatory.domains=list()

for(sp in c("human", "mouse")){
  regulatory.domains[[sp]]=list()

  datasets=system(paste("ls ", pathResults, sp, sep=""), intern=T)

  for(data in datasets){
    regulatory.domains[[sp]][[data]]=list()

    domains=system(paste("ls ", pathResults, sp, "/", data, sep=""), intern=T)

    for(domain in domains){
      regulatory.domains[[sp]][[data]][[domain]]=list()

      methods=system(paste("ls ", pathResults, sp, "/", data, "/", domain, sep=""), intern=T)
      methods=grep("GREAT", methods, value=T)

      for(method in methods){
        res=read.table(paste(pathResults, sp, "/", data, "/", domain, "/", method, "/GOntact_regulatory_domains.txt",sep=""), h=F, sep="\t", quote="\"")
        colnames(res)=c("GeneID", "Chr", "Start", "End", "Strand")

        regulatory.domains[[sp]][[data]][[domain]][[method]]=res
      }
    }
  }
}

#####################################################################################

save(regulatory.domains, file=paste(pathFigures, "RData/data.regulatory.domains.RData",sep=""))

#####################################################################################
