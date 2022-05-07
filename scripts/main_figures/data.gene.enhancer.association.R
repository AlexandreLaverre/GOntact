#####################################################################################

pathResults="../../results/"
pathFigures="../../data_for_publication/figures/"

options(stringsAsFactors=F)

#####################################################################################

for(sp in c("human", "mouse")){
  datasets=system(paste("ls ", pathResults, sp,sep=""), intern=T)

  for(data in datasets){

    domains=system(paste("ls ", pathResults, sp, "/", data, sep=""), intern=T)

    for(domain in domains){

      background.association=list()
      foreground.association=list()

      methods=system(paste("ls ", pathResults, sp, "/", data, "/", domain, sep=""), intern=T)

      for(method in methods){
        path.bg=paste(pathResults, sp, "/", data, "/", domain, "/", method, "/GOntact_element_gene_association_background.txt",sep="")

        if(file.exists(path.bg)){
          res.bg=read.table(path.bg, h=T, sep="\t", quote="\"")
        } else{
          path.bg.gz=paste(pathResults, sp, "/", data, "/", domain, "/", method, "/GOntact_element_gene_association_background.txt.gz",sep="")
          conn=gzfile(path.bg.gz, open="r")
          res.bg=read.table(conn, h=T, sep="\t", quote="\"")
          close(conn)
        }

        background.association[[method]]=res.bg

        res.fg=read.table(paste(pathResults, sp, "/", data, "/", domain, "/", method, "/GOntact_element_gene_association_foreground.txt",sep=""), h=T, sep="\t", quote="\"")

        foreground.association[[method]]=res.fg

        save(background.association, file=paste(pathFigures, "RData/data.background.association.",sp,".",data, ".",domain,".RData",sep=""))
        save(foreground.association, file=paste(pathFigures, "RData/data.foreground.association.",sp,".",data, ".", domain, ".RData",sep=""))
      }
    }
  }
}

#####################################################################################
