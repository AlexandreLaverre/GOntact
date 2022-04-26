########################################################################

pathResults="../../results/"
pathGO="../../data/GeneOntology/"
pathAnnot="../../data/ensembl_annotations/"

########################################################################

## given a list of genes, return all possible pairs

construct.pairs<-function(genes){
  genes=sort(genes)
  pairs=c()
  if(length(genes)>=2){
    for(i in 1:(length(genes)-1)){
     pairs=c(pairs, paste(genes[i], genes[(i+1):length(genes)]))
   }
  }
  return(pairs)
}

########################################################################

sp="human"
genome="hg38"
ens=102
sample="Vista_heart_vs_ENCODE"
domain="biological_process"
method.great="GREAT_upstream5kb_downstream1kb_extend1Mb"
method.contacts="contacts_mindist0kb_maxdist1Mb"

########################################################################

## gene ontology annotations

go.annot=read.table(paste(pathGO, "goa_human.gaf",sep=""), h=F, stringsAsFactors=F, sep="\t", comment.char="!", quote="\"")
go.cat=read.table(paste(pathGO, "GOCategories.txt",sep=""), h=T, stringsAsFactors=F, sep="\t", quote="\"")

this.domain=go.cat$ID[which(go.cat$GOSpace==domain)]
go.annot=go.annot[which(go.annot$V5%in%this.domain),]

########################################################################

## major tss coordinates

major.iso=read.table(paste(pathResults, sp, "/", sample, "/", domain, "/", method.great, "/GOntact_major_isoforms.txt", sep=""), h=T, stringsAsFactors=F)
annot=read.table(paste(pathAnnot,sp, "/GeneAnnotation_BioMart_Ensembl",ens,"_",genome,".txt",sep=""), h=T, stringsAsFactors=F, sep="\t", quote="\"")
annot=annot[which(annot$Transcript.stable.ID%in%majoriso$Transcript),]
rownames(annot)=annot$Transcript.stable.ID

major.iso$Chr=annot[major.iso$Transcript, "Chromosome.scaffold.name"]
major.iso$TSS=annot[major.iso$Transcript, "Transcription.start.site..TSS."]
rownames(major.iso)=major.iso$GeneID

## select only genes with GO

major.iso=major.iso[which(major.iso$GeneID%in%go.annot$V3),]

########################################################################

## bait annotations

bait.annot=read.table(paste(pathResults, sp, "/",sample, "/",domain, "/",method.contacts, "/GOntact_bait_gene_annotation.txt",sep=""), h=T, stringsAsFactors=F, sep="\t")
colnames(bait.annot)=c("BaitID", "GeneID")

## take only baits with annotated genes
bait.annot=bait.annot[which(bait.annot$GeneID!=""),]

## take only genes represented by a single bait, for simplicity
all.genes=unlist(lapply(bait.annot$GeneID, function(x) unlist(strsplit(x, split=","))))
nb.baits=as.numeric(table(as.factor(all.genes)))
names(nb.baits)=levels(as.factor(all.genes))

genes.single.bait=names(nb.baits)[which(nb.baits==1)]

## bait-gene association
nb.genes=unlist(lapply(bait.annot$GeneID, function(x) length(unlist(strsplit(x, split=",")))))

all.baits=rep(bait.annot$BaitID, nb.genes)

bait.info=data.frame("BaitID"=all.baits, "GeneID"=all.genes, stringsAsFactors=F)

########################################################################

## enhancer-gene association, contacts

enh.gene.contacts=read.table(paste(pathResults, sp, "/", sample, "/", domain, "/", method.contacts, "/GOntact_element_gene_association_background.txt", sep=""), h=T, stringsAsFactors=F, sep="\t")

## here take only genes with single bait, no ambiguity

enh.gene.contacts=enh.gene.contacts[which(enh.gene.contacts$GeneSymbol%in%genes.single.bait),]

genes.per.enh.contacts=tapply(enh.gene.contacts$GeneSymbol, as.factor(enh.gene.contacts$ElementID), function(x) x)
nb.genes.per.enh.contacts=unlist(lapply(genes.per.enh.contacts, length))

## enhancers that are in contact with multiple genes

multi.genes.contacts=genes.per.enh.contacts[which(nb.genes.per.enh.contacts>=2)]

gene.pairs.contacts=unique(unlist(lapply(multi.genes.contacts, construct.pairs)))

########################################################################

## enhancer-gene association, GREAT

enh.gene.great=read.table(paste(pathResults, sp, "/", sample, "/", domain, "/", method.great, "/GOntact_element_gene_association_background.txt", sep=""), h=T, stringsAsFactors=F, sep="\t")

genes.per.enh.great=tapply(enh.gene.great$GeneSymbol, as.factor(enh.gene.great$ElementID), function(x) x)
nb.genes.per.enh.great=unlist(lapply(genes.per.enh.great, length))

multi.genes.great=genes.per.enh.great[which(nb.genes.per.enh.great>=2)]

gene.pairs.great=unique(unlist(lapply(multi.genes.great, construct.pairs)))

########################################################################

## nb common GO annotations, contacts

contacts.nb.go.any=c()
contacts.nb.go.both=c()
contacts.id.pairs=c()
contacts.dist=c()

for(p in gene.pairs.contacts){
  s=unlist(strsplit(p, split=" "))
  g1=s[1]
  g2=s[2]
  b1=bait.info$BaitID[which(bait.info$GeneID==g1)]
  b2=bait.info$BaitID[which(bait.info$GeneID==g2)]

  if(b1!=b2){
    contacts.id.pairs=c(contacts.id.pairs, p)
    contacts.dist=c(contacts.dist, abs(major.iso[g1,"TSS"]-major.iso[g2,"TSS"]))
    go1=go.annot$V5[which(go.annot$V3==g1)]
    go2=go.annot$V5[which(go.annot$V3==g2)]
    contacts.nb.go.any=c(contacts.nb.go.any, length(unique(c(go1,go2))))
    contacts.nb.go.both=c(contacts.nb.go.both, length(intersect(go1,go2)))

    if(length(contacts.id.pairs)%%1000==0){
      print(length(contacts.id.pairs))
    }
  }
}

contacts.class.dist=cut(contacts.dist, breaks=seq(from=0, to=1e6, by=50e3), include.lowest=T)

########################################################################

## nb common GO annotations, great

great.nb.go.any=c()
great.nb.go.both=c()
great.id.pairs=c()
great.dist=c()

for(p in gene.pairs.great){
  s=unlist(strsplit(p, split=" "))
  g1=s[1]
  g2=s[2]

  great.dist=c(great.dist, abs(major.iso[g1,"TSS"]-major.iso[g2,"TSS"]))
  great.id.pairs=c(great.id.pairs, p)
  go1=go.annot$V5[which(go.annot$V3==g1)]
  go2=go.annot$V5[which(go.annot$V3==g2)]
  great.nb.go.any=c(great.nb.go.any, length(unique(c(go1,go2))))
  great.nb.go.both=c(great.nb.go.both, length(intersect(go1,go2)))

  if(length(great.id.pairs)%%1000==0){
    print(length(great.id.pairs))
  }
}

 great.class.dist=cut(great.dist, breaks=seq(from=0, to=1e6, by=50e3), include.lowest=T)

########################################################################

## nb common GO annotations, random

random.nb.go.any=c()
random.nb.go.both=c()
random.id.pairs=c()
random.dist=c()

for(i in 1:50000){
  g1=sample(major.iso$GeneID, size=1)
  chr=major.iso[g1, "Chr"]
  g2=sample(setdiff(major.iso$GeneID[which(major.iso$Chr==chr)], g1), size=1)

  random.dist=c(random.dist, abs(major.iso[g1,"TSS"]-major.iso[g2,"TSS"]))
  random.id.pairs=c(random.id.pairs, paste(g1, g2))
  go1=go.annot$V5[which(go.annot$V3==g1)]
  go2=go.annot$V5[which(go.annot$V3==g2)]
  random.nb.go.any=c(random.nb.go.any, length(unique(c(go1,go2))))
  random.nb.go.both=c(random.nb.go.both, length(intersect(go1,go2)))

  if(length(random.id.pairs)%%1000==0){
    print(length(random.id.pairs))
  }
}

random.class.dist=cut(random.dist, breaks=seq(from=0, to=1e6, by=50e3), include.lowest=T)

########################################################################
