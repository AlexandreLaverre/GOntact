#############################################################################

pathFigures="../../data_for_publication/figures/"
pathResults="../../results/motif_enrichment/"

#############################################################################

sp="mouse"
datasets=c("HNF4a", "Hoxd13", "Nr2f1_up")
motifs=c("MA0114.3", "MA0909.3", "MA1537.1")

#############################################################################

bg.enh=read.table(paste(pathResults, sp, "/ENCODE.Laverre2022.random50K.bed",sep=""), h=F, sep="\t")

nb.with=list()
nb.tot=list()
scores=list()

for(i in 1:length(motifs)){
  motif=motifs[i]
  dataset=datasets[i]

  bg.fimo=read.table(paste(pathResults, sp, "/", motif, "/fimo.tsv",sep=""), h=T, stringsAsFactors=F)

  nb.with[[motif]]=c(length(unique(bg.fimo[,3])))
  nb.tot[[motif]]=c(nrow(bg.enh))

  scores[[motif]][["bg"]]=bg.fimo$score

  for(type in c("GREAT_only", "GOntact_only", "shared")){
    fg.enh=read.table(paste(pathResults, sp, "/", dataset, "/", type, "/associated_enhancers.bed",sep=""), h=F, stringsAsFactors=F)
    fg.fimo=read.table(paste(pathResults, sp, "/", dataset, "/", type, "/",motif,"/fimo.tsv",sep=""), h=T, stringsAsFactors=F)

    nb.tot[[motif]]=c(nb.tot[[motif]], nrow(fg.enh))
    nb.with[[motif]]=c(nb.with[[motif]], length(unique(fg.fimo[,3])))

    scores[[motif]][[type]]=fg.fimo$score

  }

  names(nb.tot[[motif]])=c("background", "GREAT_only", "GOntact_only", "shared")
  names(nb.with[[motif]])=c("background", "GREAT_only", "GOntact_only", "shared")
}

#############################################################################

save(list=c("nb.tot", "nb.with", "scores"), file=paste(pathFigures, "RData/data.motif.enrichment.RData",sep=""))

#############################################################################
