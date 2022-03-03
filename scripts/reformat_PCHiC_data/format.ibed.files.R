###########################################################

path="../../data/PCHi-C/"

###########################################################

for(sp in c("human", "mouse")){

  files=system(paste("ls ",path,sp, "/processed_ibed_files/ | grep ibed",sep=""), intern=T)

  samples=unlist(lapply(files, function(x) unlist(strsplit(x,split="\\."))[1]))

  for(sample in samples){
    print(sample)
    coords=read.table(paste(path, sp, "/processed_ibed_files/",sample,".ibed",sep=""),h=T, stringsAsFactors=F, sep="\t")

    coords$bait_chr=unlist(lapply(coords$bait_chr, function(x) substr(x,4, nchar(x))))
    coords$chr=unlist(lapply(coords$chr, function(x) substr(x,4, nchar(x))))

    ## only standard chromosomes
    coords=coords[which(coords$bait_chr%in%c(as.character(1:22), "X", "Y") & coords$chr%in%c(as.character(1:22), "X", "Y")),]

    coords$bait_name=paste(coords$bait_chr, ":",coords$bait_start, "-", coords$bait_end,sep="")
    coords$otherEnd_name=paste(coords$chr, ":",coords$start, "-", coords$end,sep="")

    coords=coords[,c("bait_chr", "bait_start", "bait_end", "bait_name", "chr", "start", "end", "otherEnd_name", "N_reads", "score")]
    colnames(coords)=c("bait_chr", "bait_start", "bait_end", "bait_name", "otherEnd_chr", "otherEnd_start", "otherEnd_end", "otherEnd_name", "N_reads", "score")

    write.table(coords, file=paste(path, sp, "/original_ibed_files/",sample,".ibed",sep=""), row.names=F, col.names=T, sep="\t", quote=F)

  }
}

###########################################################
