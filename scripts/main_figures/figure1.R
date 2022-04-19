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
  all.enhancers$ID=paste(all.enhancers$Chr,":", all.enhancers$Start,"-",all.enhancers$End,sep="")
  all.enhancers$Great=all.enhancers$ID%in%enhancers.great$ElementID
  all.enhancers$Contacts=all.enhancers$ID%in%enhancers.contacts$ElementID
  
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

for(i in 1:5){
  m[i,]=rep(1,5)
}

for(i in 5:6){
  m[i,]=rep(2,5)
}

for(i in 7:12){
  m[i,]=rep(3,5)
}

layout(m)

#####################################################################################

## plot gene major TSS

ylim=c(-0.4,1)

par(mar=c(0.1, 2.1, 2.1, 1))

plot(1, type="n", xlab="", ylab="", axes=F, main="", xaxs="i", yaxs="i", xlim=xlim, ylim=ylim)

segments(xlim[1], 0, xlim[2], 0)

gene.colors=rep("black",length(all.tss))
names(gene.colors)=names(all.tss)
gene.colors["POU6F1"]="red"

tinyy=0.02

ypos.domains=c(-0.45, -0.65, -0.9, -0.45, -0.65, -0.9)
names(ypos.domains)=names(all.tss)

for(g in names(all.tss)){
  this.tss=all.tss[g]
  this.strand=all.strands[g]
  this.col=gene.colors[g]

  if(this.strand=="+"){
    height=0.15
    arrowsize=diff(xlim)/25

    basal.start=this.tss-5000
    basal.end=this.tss+1000
  }

  if(this.strand=="-"){
    height=-0.15
    arrowsize=-diff(xlim)/25
    basal.start=this.tss+5000
    basal.end=this.tss-1000
  }

  rect(basal.start, -tinyy, basal.end, tinyy, col=this.col, border=this.col)

  segments(this.tss, 0, this.tss, height, col=this.col)
  arrows(this.tss, height, this.tss+arrowsize, height, length=0.05, col=this.col)

  if(g=="POU6F1"){
    text(g, x=this.tss-diff(xlim)/30, y=0.05, adj=c(0.025,0.02), font=3, xpd=NA, col=this.col)
   }
}

text("TSS", x=xlim[2]-diff(xlim)/40, y=0.1, cex=1.1)

#####################################################################################

## draw enhancers on the same plot

yenh=0.75
height.enh=0.025
segments(xlim[1], yenh, xlim[2], yenh)

text("enhancers", x=xlim[2]-diff(xlim)/22, y=0.9, cex=1.1)

segments(all.enhancers$Start[which(!(all.enhancers$Great | all.enhancers$Contacts))] , yenh-height.enh, all.enhancers$End[which(!(all.enhancers$Great | all.enhancers$Contacts))], yenh+height.enh)
segments(all.enhancers$Start[which(all.enhancers$Great)] , yenh-height.enh, all.enhancers$End[which(all.enhancers$Great)], yenh+height.enh, col="red")
segments(all.enhancers$Start[which(all.enhancers$Contacts)] , yenh-height.enh, all.enhancers$End[which(all.enhancers$Contacts)], yenh+height.enh, col="steelblue")

#####################################################################################

## regulatory domains on the same plot

gene.ypos=c(0.55, 0.41, 0.3, 0.55, 0.41, 0.3)
names(gene.ypos)=names(all.tss)

text("domains", x=xlim[2]-diff(xlim)/28, y=0.38, cex=1.1)

for(g in names(all.tss)){
  this.domain.start=this.domains$Start[which(this.domains$GeneID==g)]
  this.domain.end=this.domains$End[which(this.domains$GeneID==g)]

  this.col=gene.colors[g]

  segments(this.domain.start, gene.ypos[g], this.domain.end, gene.ypos[g], col=this.col)
  
  if(g=="POU6F1"){
    segments(this.domain.start, 0, this.domain.start, 0.75, lty=3, col="red")
    segments(this.domain.end, 0, this.domain.end, 0.75, lty=3, col="red")
  }
}

prettyx=pretty(xlim)
labels=paste(prettyx/1e6, "Mb", sep=" ")
axis(side=1, at=prettyx, labels=labels, mgp=c(3,0.5,0))

#####################################################################################

## contacted enhancers

for(i in which(all.enhancers$Contacts)){
  pos.enh=(all.enhancers$Start[i]+all.enhancers$End[i])/2
  pos.gene=all.tss["POU6F1"]
  x.center=(pos.enh+pos.gene)/2

  y.first=0.82
  y.top=1.15
  
  segments(pos.gene, y.first, x.center, y.top, col="lightblue", xpd=NA)
  segments(x.center, y.top, pos.enh, y.first, col="lightblue", xpd=NA)
}

## plot label

mtext("A", side=3, line=0.1, at=xlim[1]-diff(xlim)/50, font=2, cex=1.1)

#####################################################################################

dev.off()

#####################################################################################
