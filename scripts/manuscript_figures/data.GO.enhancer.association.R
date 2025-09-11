#####################################################################################

pathResults="../../results/"
pathFigures="../../data_for_publication/figures/"

options(stringsAsFactors=F)

#####################################################################################

foreground.association=list()

for(sp in c("human", "mouse")){
  foreground.association[[sp]]=list()

  datasets=system(paste("ls ", pathResults, sp,sep=""), intern=T)

  for(data in datasets){
    foreground.association[[sp]][[data]]=list()

    domains=system(paste("ls ", pathResults, sp, "/", data, sep=""), intern=T)

    for(domain in domains){
      foreground.association[[sp]][[data]][[domain]]=list()

      methods=system(paste("ls ", pathResults, sp, "/", data, "/", domain," | grep -v txt", sep=""), intern=T)

      for(method in methods){
        res.fg=read.table(paste(pathResults, sp, "/", data, "/", domain, "/", method, "/GOntact_element_GO_association_foreground.txt",sep=""), h=F, sep="\t", quote="\"")

        foreground.association[[sp]][[data]][[domain]][[method]]=res.fg
      }
    }
  }
}

#####################################################################################

save(foreground.association, file=paste(pathFigures, "RData/data.foreground.GO.association.RData",sep=""))

#####################################################################################
