#####################################################################################

pathFigures="../../data_for_publication/figures/"

options(stringsAsFactors=F)

library(plotrix)

#####################################################################################

if(load){
  load(paste(pathFigures, "RData/data.background.association.RData",sep=""))
  load(paste(pathFigures, "RData/data.regulatory.domains.RData",sep=""))

  load(paste(pathFigures, "RData/data.ENCODE.human.RData",sep=""))
  colnames(encode)=c("Chr", "Start", "End")

  sp="human"
  dataset="Vista_heart_vs_ENCODE"
  domain="biological_process"

  load=FALSE
}

#####################################################################################

if(prepare){

  assoc.great=background.association[[sp]][[dataset]][[domain]][["GREAT_upstream5kb_downstream1kb_extend1Mb"]]
  assoc.contacts=background.association[[sp]][[dataset]][[domain]][["contacts_mindist25kb_maxdist1Mb"]]

  ## enhancer coordinates, GREAT association
  assoc.great$Chr=unlist(lapply(assoc.great$ElementID, function(x) unlist(strsplit(x, split=":"))[1]))
  assoc.great$Coords=unlist(lapply(assoc.great$ElementID, function(x) unlist(strsplit(x, split=":"))[2]))
  assoc.great$Start=as.numeric(unlist(lapply(assoc.great$Coords, function(x) unlist(strsplit(x, split="-"))[1])))
  assoc.great$End=as.numeric(unlist(lapply(assoc.great$Coords, function(x) unlist(strsplit(x, split="-"))[2])))

  ## enhancer coordinates, contacts

  assoc.contacts$Chr=unlist(lapply(assoc.contacts$ElementID, function(x) unlist(strsplit(x, split=":"))[1]))
  assoc.contacts$Coords=unlist(lapply(assoc.contacts$ElementID, function(x) unlist(strsplit(x, split=":"))[2]))
  assoc.contacts$Start=as.numeric(unlist(lapply(assoc.contacts$Coords, function(x) unlist(strsplit(x, split="-"))[1])))
  assoc.contacts$End=as.numeric(unlist(lapply(assoc.contacts$Coords, function(x) unlist(strsplit(x, split="-"))[2])))

  ## regulatory domains
  domains=regulatory.domains[[sp]][[dataset]][[domain]][["GREAT_upstream5kb_downstream1kb_extend1Mb"]]
  domains$Size=domains$End-domains$Start+1

  ## take only genes that have at least 1 enhancer in both datasets
  common.genes=intersect(assoc.great$GeneSymbol, assoc.contacts$GeneSymbol)

  assoc.great.common=assoc.great[which(assoc.great$GeneSymbol%in%common.genes),]
  assoc.contacts.common=assoc.contacts[which(assoc.contacts$GeneSymbol%in%common.genes),]

  nbenh.great=table(assoc.great.common$GeneSymbol)
  nbenh.contacts=table(assoc.contacts.common$GeneSymbol)

  mediandist.great=tapply(assoc.great.common$Distance, as.factor(assoc.great.common$GeneSymbol), median)
  mediandist.contacts=tapply(assoc.contacts.common$Distance, as.factor(assoc.contacts.common$GeneSymbol), median)

  ## example gene

  gene="POU6F1"
  this.tss=51217671
  this.start.domain=domains$Start[which(domains$GeneID==gene)]
  this.end.domain=domains$End[which(domains$GeneID==gene)]
  this.chr=domains$Chr[which(domains$GeneID==gene)]
  this.strand="-"

  enhancers.great=assoc.great.common[which(assoc.great.common$GeneSymbol==gene),]
  enhancers.contacts=assoc.contacts.common[which(assoc.contacts.common$GeneSymbol==gene),]

  ## axis limits

  xlim=c(min(c(enhancers.great$Start, enhancers.contacts$Start)), max(c(enhancers.great$End, enhancers.contacts$End)))

  this.domains=domains[which(domains$Chr==this.chr & ((domains$Start>=xlim[1] & domains$Start<=xlim[2]) | (domains$End>=xlim[1] & domains$End<=xlim[2]))),]
  xlim[1]=min(c(xlim[1], this.domains$Start))
  xlim[2]=max(c(xlim[2], this.domains$End))

  all.tss=c(50953922, 51026384, 51083596, 51173135, 51217671, 51238826)
  names(all.tss)=c("HIGD1C", "SLC11A2", "CSRNP2", "TFCP2", "POU6F1", "DAZAP2")

  all.strands=c("+", "-", "-", "-", "-", "+")
  names(all.strands)=c("HIGD1C", "SLC11A2", "CSRNP2", "TFCP2", "POU6F1", "DAZAP2")

  all.enhancers=encode[which(encode$Chr==paste("chr",this.chr,sep="") & encode$Start>=xlim[1] & encode$End<=xlim[2]),]

  prepare=FALSE
}

#####################################################################################
#####################################################################################

## 1 column width 85 mm = 3.34 in
## 1.5 column width 114 mm = 4.49 in
## 2 columns width 174 mm = 6.85 in
## max height: 11 in

pdf(file=paste(pathFigures, "Figure1.pdf",sep=""), width=6.85, height=5.5)

m=matrix(rep(NA, 12*5), nrow=12)

for(i in 1:3){
  m[i,]=rep(1,5)
}

for(i in 4:5){
  m[i,]=rep(2,5)
}

for(i in 6:12){
  m[i,]=rep(3,5)
}

layout(m)

#####################################################################################

## plot gene major TSS

ylim=c(-0.3,0.45)

par(mar=c(1,1,2.1,1))

plot(1, type="n", xlab="", ylab="", axes=F, main="", xaxs="i", yaxs="i", xlim=xlim, ylim=ylim)

segments(xlim[1], 0, xlim[2], 0)

gene.colors=rep("black",length(all.tss))
names(gene.colors)=names(all.tss)
gene.colors["POU6F1"]="red"

tinyy=0.02

ypos.domains=c(-0.485, -0.72, -0.95, -0.485, -0.72, -0.95)
names(ypos.domains)=names(all.tss)

for(g in names(all.tss)){
  this.tss=all.tss[g]
  this.strand=all.strands[g]
  this.col=gene.colors[g]

  if(this.strand=="+"){
    height=0.25
    arrowsize=diff(xlim)/25

    basal.start=this.tss-5000
    basal.end=this.tss+1000
  }

  if(this.strand=="-"){
    height=-0.25
    arrowsize=-diff(xlim)/25
    basal.start=this.tss+5000
    basal.end=this.tss-1000
  }

  rect(basal.start, -tinyy, basal.end, tinyy, col=this.col, border=this.col)

  segments(this.tss, 0, this.tss, height, col=this.col)
  arrows(this.tss, height, this.tss+arrowsize, height, length=0.05, col=this.col)

  segments(this.tss, 0, this.tss,ypos.domains[g], col=this.col, xpd=NA, lty=3)
}

prettyx=pretty(xlim)
labels=paste(prettyx/1e6, "Mb", sep=" ")
axis(side=3, at=prettyx, labels=labels, mgp=c(3,0.5,0))

#####################################################################################

## plot regulatory domains

## gene.ypos=1:length(all.tss)

gene.ypos=c(6,5,4,6,5,4)
names(gene.ypos)=names(all.tss)

par(mar=c(1,1,0.25,1))

plot(1, type="n", xlab="", ylab="", axes=F, main="", xaxs="i", yaxs="i", xlim=xlim, ylim=c(3.75,6.25))

for(g in names(all.tss)){
  this.domain.start=this.domains$Start[which(this.domains$GeneID==g)]
  this.domain.end=this.domains$End[which(this.domains$GeneID==g)]

  this.col=gene.colors[g]

  segments(this.domain.start, gene.ypos[g], this.domain.end, gene.ypos[g], col=this.col)

  text(g, x=this.domain.start, y=gene.ypos[g], adj=c(0,-0.2), font=3, xpd=NA)
}

#####################################################################################

#####################################################################################

dev.off()

#####################################################################################
