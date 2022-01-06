##########################################################################

pathEnsembl="../../data/ensembl_annotations/"

release=94

##########################################################################

library(seqinr)

##########################################################################

for(sp in c("human", "mouse")){
  ## read transcript CDS sequences

  cds=read.fasta(paste(pathEnsembl, sp, "/AllCDS_Ensembl", release, ".fa",sep=""), seqtype="DNA", forceDNAtolower=FALSE)

  lencds=unlist(lapply(cds, length))
  names(lencds)=unlist(lapply(names(lencds), function(x) unlist(strsplit(x, split="\\."))[1]))

  ## APPRIS annotations
  appris=read.table(paste(pathEnsembl, sp, "/APPRIS_Ensembl",release,".txt", sep=""), h=T, stringsAsFactors=F, sep="\t")
  rownames(appris)=appris$Transcript.stable.ID

  ## transcript info
  txinfo=read.table(paste(pathEnsembl, sp, "/TranscriptInfo_Ensembl", release,".txt", sep=""), h=T, stringsAsFactors=F, sep="\t", quote="\"")
  colnames(txinfo)[1]="gene_id"
  colnames(txinfo)[2]="transcript_id"
  colnames(txinfo)[4]="chr"
  colnames(txinfo)[5]="start"
  colnames(txinfo)[6]="end"
  colnames(txinfo)[7]="strand"
  
  ## gene info
  geneinfo=read.table(paste(pathEnsembl, sp, "/GeneInfo_Ensembl", release,".txt", sep=""), h=T, stringsAsFactors=F, sep="\t", quote="\"")
  rownames(geneinfo)=geneinfo$stable_id

  ## add gene biotype to transcript info
  txinfo$gene_biotype=geneinfo[txinfo$gene_id,"biotype"]

  ## we take only protein-coding genes
  ## select protein-coding transcripts
   
  txinfo=txinfo[which(txinfo$biotype=="protein_coding" & txinfo$gene_biotype=="protein_coding"),]

  ## add transcript CDS length
  
  txinfo$length=lencds[txinfo$transcript_id]
  
  ## remove NA transcript lengths - haplotypes

  txinfo=txinfo[which(!is.na(txinfo$length)),]

  ## order transcripts by length

  txinfo=txinfo[order(txinfo$length, decreasing=T),]

  ## put "principal transcripts" first

  txinfo$APPRIS=appris[txinfo$transcript_id, "APPRIS.annotation"]
  p1=which(txinfo$APPRIS=="principal1")
  p2=which(txinfo$APPRIS=="principal2")
  p3=which(txinfo$APPRIS=="principal3")
  p4=which(txinfo$APPRIS=="principal4")
  p5=which(txinfo$APPRIS=="principal5")
  a1=which(txinfo$APPRIS=="alternative1")
  a2=which(txinfo$APPRIS=="alternative2")
  
  other=setdiff(1:dim(txinfo)[1], c(p1, p2, p3, p4, p5, a1, a2))
  txinfo=txinfo[c(p1, p2, p3, p4, p5, a1, a2, other),]

  ## select the first transcript for each gene
  ## APPRIS principal for the genes that have them
  ## longest CDS for the other ones

  first=which(!duplicated(txinfo$gene_id))
  txinfo=txinfo[first,]

  ## extract TSS

  txinfo$TSS=rep(NA, dim(txinfo)[1])
  txinfo$TSS[which(txinfo$strand==1)]=txinfo$start[which(txinfo$strand==1)]
  txinfo$TSS[which(txinfo$strand==-1)]=txinfo$end[which(txinfo$strand==-1)]

  write.table(txinfo, file=paste(pathEnsembl, sp, "/canonical_transcripts_Ensembl",release,"_protein_coding.txt",sep=""), row.names=F, col.names=T, sep="\t", quote=F)
  
}

##########################################################################
