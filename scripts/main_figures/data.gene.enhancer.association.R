#####################################################################################

pathResults="../../results/"
pathFigures="../../data_for_publication/figures/"

options(stringsAsFactors=F)

#####################################################################################

background.association=list()
foreground.association=list()

for(sp in c("human", "mouse")){
  background.association[[sp]]=list()
  foreground.association[[sp]]=list()

  datasets=system(paste("ls ", pathResults, sp,sep=""), intern=T)

  for(data in datasets){
    background.association[[sp]][[data]]=list()
    foreground.association[[sp]][[data]]=list()

    domains=system(paste("ls ", pathResults, sp, "/", data, sep=""), intern=T)

    for(domain in domains){
      background.association[[sp]][[data]][[domain]]=list()
      foreground.association[[sp]][[data]][[domain]]=list()

      methods=system(paste("ls ", pathResults, sp, "/", data, "/", domain, sep=""), intern=T)

      for(method in methods){
        res.bg=read.table(paste(pathResults, sp, "/", data, "/", domain, "/", method, "/GOntact_element_gene_association_background.txt",sep=""), h=T, sep="\t", quote="\"")

        background.association[[sp]][[data]][[domain]][[method]]=res.bg

        res.fg=read.table(paste(pathResults, sp, "/", data, "/", domain, "/", method, "/GOntact_element_gene_association_foreground.txt",sep=""), h=T, sep="\t", quote="\"")

        foreground.association[[sp]][[data]][[domain]][[method]]=res.fg
      }
    }
  }
}

#####################################################################################

save(background.association, file=paste(pathFigures, "RData/data.background.association.RData",sep=""))
save(foreground.association, file=paste(pathFigures, "RData/data.foreground.association.RData",sep=""))

#####################################################################################
