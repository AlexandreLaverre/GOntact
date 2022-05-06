###########################################################

path="../../data/PCHi-C/"

###########################################################

for(sp in c("human", "mouse")){

  files=system(paste("ls ",path,sp, "/simulated_data/ | grep txt",sep=""), intern=T)

  samples=unlist(lapply(files, function(x) unlist(strsplit(x,split="\\."))[1]))
  samples=unlist(lapply(samples, function(x) paste(unlist(strsplit(x, split="_"))[-1], collapse="_")))

  for(sample in samples){
    print(sample)
    coords=read.table(paste(path, sp, "/simulated_data/simulated_",sample,".txt",sep=""),h=T, stringsAsFactors=F, sep="\t")
    coords=coords[,1:8]

    coords$bait_chr=unlist(lapply(coords$bait_chr, function(x) substr(x,4, nchar(x))))
    coords$chr=unlist(lapply(coords$chr, function(x) substr(x,4, nchar(x))))

    ## only standard chromosomes
    coords=coords[which(coords$bait_chr%in%c(as.character(1:22), "X", "Y") & coords$chr%in%c(as.character(1:22), "X", "Y")),]

    coords$bait_name=paste(coords$bait_chr, ":",coords$bait_start, "-", coords$bait_end,sep="")
    coords$otherEnd_name=paste(coords$chr, ":",coords$start, "-", coords$end,sep="")

    coords=coords[,c("bait_chr", "bait_start", "bait_end", "bait_name", "chr", "start", "end", "otherEnd_name", "N_reads", "score")]
    colnames(coords)=c("bait_chr", "bait_start", "bait_end", "bait_name", "otherEnd_chr", "otherEnd_start", "otherEnd_end", "otherEnd_name", "N_reads", "score")

    coords$N_reads=rep(0, nrow(coords))
    coords$score=rep(0, nrow(coords))

    write.table(coords, file=paste(path, sp, "/simulated_ibed_files/",sample,".ibed",sep=""), row.names=F, col.names=T, sep="\t", quote=F)

  }
}

###########################################################
