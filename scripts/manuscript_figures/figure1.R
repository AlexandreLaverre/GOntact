#####################################################################################

pathFigures="../../data_for_publication/figures/"

options(stringsAsFactors=F)

#####################################################################################

if(load){

  load(paste(pathFigures, "RData/data.regulatory.domains.RData",sep=""))
  load(paste(pathFigures, "RData/data.ENCODE.human.RData",sep=""))
  colnames(encode)=c("Chr", "Start", "End")

  sp="human"
  dataset="Vista_heart_vs_ENCODE"
  domain="biological_process"

  load(paste(pathFigures, "RData/data.background.association.",sp,".",dataset, ".",domain,".RData",sep=""))

  load=FALSE
}

#####################################################################################

if(prepare){
  ## regulatory domains
  domains=regulatory.domains[[sp]][[dataset]][[domain]][["GREAT_upstream5kb_downstream1kb_extend1Mb"]]
  domains$Size=domains$End-domains$Start+1

  ## gene-enhancer association, all ENCODE

  results=list()

  methods=c("GREAT_upstream5kb_downstream1kb_extend1Mb", "GREAT_fixed_size_upstream1Mb_downstream1Mb", "contacts_mindist0kb_maxdist1Mb",  "shared_contacts_mindist0kb_maxdist1Mb", "hybrid_mindist25kb_maxdist1Mb")
  shortnames=c("GREAT", "fixed size", "GOntact, all data", "GOntact, shared contacts", "GOntact, hybrid")

  names(shortnames)=methods

  for(method in methods){
    print(method)

    assoc=background.association[[method]]

    ## enhancer coordinates

    if(method=="GREAT_upstream5kb_downstream1kb_extend1Mb" | method=="contacts_mindist0kb_maxdist1Mb"){
      assoc$Chr=unlist(lapply(assoc$ElementID, function(x) unlist(strsplit(x, split=":"))[1]))
      assoc$Coords=unlist(lapply(assoc$ElementID, function(x) unlist(strsplit(x, split=":"))[2]))
      assoc$Start=as.numeric(unlist(lapply(assoc$Coords, function(x) unlist(strsplit(x, split="-"))[1])))
      assoc$End=as.numeric(unlist(lapply(assoc$Coords, function(x) unlist(strsplit(x, split="-"))[2])))
    }

    nb.enh=table(assoc$GeneSymbol)

    nb.genes=table(assoc$ElementID)

    assoc$DistanceClass=cut(assoc$Distance, breaks=seq(from=0, to=1e6, by=50e3), include.lowest=T)

    prop.dist=as.numeric(table(assoc$DistanceClass))/sum(as.numeric(table(assoc$DistanceClass)))

    if(method=="GREAT_upstream5kb_downstream1kb_extend1Mb" | method=="contacts_mindist0kb_maxdist1Mb"){
      results[[method]]=list("assoc"=assoc, "nb.enh"=nb.enh, "nb.genes"=nb.genes, "prop.dist"=prop.dist)
    } else{
      results[[method]]=list("nb.enh"=nb.enh, "nb.genes"=nb.genes, "prop.dist"=prop.dist)
    }
  }

  assoc.great=results[["GREAT_upstream5kb_downstream1kb_extend1Mb"]][["assoc"]]
  assoc.contacts=results[["contacts_mindist0kb_maxdist1Mb"]][["assoc"]]

  ## example gene

  gene="POU6F1"
  this.tss=51217671
  this.start.domain=domains$Start[which(domains$GeneID==gene)]
  this.end.domain=domains$End[which(domains$GeneID==gene)]
  this.chr=domains$Chr[which(domains$GeneID==gene)]
  this.strand="-"

  enhancers.great=assoc.great[which(assoc.great$GeneSymbol==gene),]
  enhancers.contacts=assoc.contacts[which(assoc.contacts$GeneSymbol==gene),]

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

m=matrix(rep(NA, 12*10), nrow=12)

for(i in 1:5){
  m[i,]=rep(1,10)
}

for(i in 6:12){
  m[i,]=c(rep(2,3), rep(3, 2),  rep(4,5))
}

layout(m)

#####################################################################################

## plot gene major TSS

ylim=c(-0.4,1)

par(mar=c(2.1, 1.5, 2.1, 1.5))

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
segments(all.enhancers$Start[which(all.enhancers$Contacts)] , yenh-height.enh, all.enhancers$End[which(all.enhancers$Contacts)], yenh+height.enh, col="darkorange")

text("PCHi-C contacts", x=xlim[2]-diff(xlim)/15, y=1.05, cex=1.1, col="darkorange", xpd=NA)

#####################################################################################

## regulatory domains on the same plot

gene.ypos=c(0.55, 0.41, 0.3, 0.55, 0.41, 0.3)
names(gene.ypos)=names(all.tss)

text("GREAT domains", x=xlim[2]-diff(xlim)/14, y=0.38, cex=1.1)

for(g in names(all.tss)){
  this.domain.start=this.domains$Start[which(this.domains$GeneID==g)]
  this.domain.end=this.domains$End[which(this.domains$GeneID==g)]

  this.col=gene.colors[g]

  segments(this.domain.start, gene.ypos[g], this.domain.end, gene.ypos[g], col=this.col)

  if(g=="POU6F1"){
    segments(this.domain.start, 0, this.domain.start, 0.75, lty=3, col="black")
    segments(this.domain.end, 0, this.domain.end, 0.75, lty=3, col="black")
  }
}

prettyx=pretty(xlim)
labels=paste(prettyx/1e6, "Mb", sep=" ")
axis(side=1, at=prettyx, labels=labels, mgp=c(3,0.5,0), cex.axis=1.05)
mtext("chr12", side=1, line=0.5, at=xlim[1]+diff(xlim)/50, cex=0.75)

#####################################################################################

## contacted enhancers

for(i in which(all.enhancers$Contacts)){
  pos.enh=(all.enhancers$Start[i]+all.enhancers$End[i])/2
  pos.gene=all.tss["POU6F1"]
  x.center=(pos.enh+pos.gene)/2

  y.first=0.82
  y.top=1.15

  segments(pos.gene, y.first, x.center, y.top, col="darkorange", xpd=NA)
  segments(x.center, y.top, pos.enh, y.first, col="darkorange", xpd=NA)
}

## plot label

mtext("A", side=3, line=0.5, at=xlim[1]-diff(xlim)/55, font=2, cex=1.1)

#####################################################################################

## boxplot for the number of enhancers per gene

par(mar=c(5.1, 3.5, 2.1, 1.5))
nb.enh.great=results[["GREAT_upstream5kb_downstream1kb_extend1Mb"]][["nb.enh"]]
nb.enh.contacts=results[["contacts_mindist0kb_maxdist1Mb"]][["nb.enh"]]
nb.enh.shared.contacts=results[["shared_contacts_mindist0kb_maxdist1Mb"]][["nb.enh"]]
nb.enh.hybrid=results[["hybrid_mindist25kb_maxdist1Mb"]][["nb.enh"]]
nb.enh.fixed=results[["GREAT_fixed_size_upstream1Mb_downstream1Mb"]][["nb.enh"]]

boxplot(as.numeric(nb.enh.great), as.numeric(nb.enh.contacts), as.numeric(nb.enh.shared.contacts), as.numeric(nb.enh.hybrid),  col="white", border=c("red", "darkorange", "slateblue", "steelblue"),  outline=F, names=rep("", 2), lwd=1.25, boxwex=0.75, notch=T, axes=F, ylim=c(0,200))
axis(side=2, mgp=c(3,0.75,0))
axis(side=1, at=1:4, labels=rep("",4))

mtext("nb. enhancers per gene", side=2, line=2.35, cex=0.75)

legend("bottomleft", legend=c("GREAT (basal + extension)", "GOntact (all PCHi-C data)", "GOntact (common contacts)"), col=c("red", "darkorange", "slateblue"), lty=1, bty="n", inset=c(-0.35,-0.3), xpd=NA)
legend("bottomleft", legend=c("GOntact (hybrid 25kb)","fixed 2Mb window"), col=c("steelblue", "black"), lty=1, bty="n", inset=c(1,-0.25), xpd=NA)

## plot label

mtext("B", side=3, line=0.5, at=-0.75, font=2, cex=1.1)

## fixed size

par(mar=c(5.1, 3.5, 2.1, 1))
boxplot(as.numeric(nb.enh.fixed),  col="white", border=c("black"),  outline=F, names=rep("", 1), lwd=1.25, boxwex=0.75, notch=T, axes=F)
axis(side=2, mgp=c(3,0.75,0))
axis(side=1, at=c(0.25,1,1.75), labels=rep("",3))
mtext("nb. enhancers per gene", side=2, line=2.35, cex=0.75)

mtext("C", side=3, line=0.5, at=-0.05, font=2, cex=1.1)

#####################################################################################

## distance distribution

prop.dist.great=results[["GREAT_upstream5kb_downstream1kb_extend1Mb"]][["prop.dist"]]
prop.dist.contacts=results[["contacts_mindist0kb_maxdist1Mb"]][["prop.dist"]]
prop.dist.shared.contacts=results[["shared_contacts_mindist0kb_maxdist1Mb"]][["prop.dist"]]
prop.dist.fixed=results[["GREAT_fixed_size_upstream1Mb_downstream1Mb"]][["prop.dist"]]
prop.dist.hybrid=results[["hybrid_mindist25kb_maxdist1Mb"]][["prop.dist"]]

xpos=1:20
this.xlim=c(0.5,20.5)
this.ylim=c(0, 40)

par(mar=c(4.1,4.1,2.1,1.1))

plot(1, type="n", xlab="", ylab="", xlim=this.xlim, ylim=this.ylim, main="", axes=F)

lines(xpos, 100*prop.dist.great, col="red")
lines(xpos, 100*prop.dist.contacts, col="darkorange")
lines(xpos, 100*prop.dist.shared.contacts, col="slateblue")
lines(xpos, 100*prop.dist.fixed,  col="black")
lines(xpos, 100*prop.dist.hybrid,  col="steelblue")

points(xpos, 100*prop.dist.great, pch=20, col="red", cex=1.1)
points(xpos, 100*prop.dist.contacts, pch=20, col="darkorange", cex=1.1)
points(xpos, 100*prop.dist.shared.contacts, pch=20, col="slateblue", cex=1.1)
points(xpos, 100*prop.dist.fixed, pch=20, col="black", cex=1.1)
points(xpos, 100*prop.dist.hybrid, pch=20, col="steelblue", cex=1.1)


xax=c(1, 5.5, 10.5, 15.5, 20)
xax.lab=c("25kb", "250kb", "500kb", "750kb", "1Mb")
axis(side=1, at=xax, labels=xax.lab)
axis(side=2, cex.axis=1.05)

mtext("% gene-enhancer associations", side=2, line=2.75, cex=0.75)
mtext("enhancer - TSS distance", side=1, line=2.5, cex=0.75)

## plot label

mtext("D", side=3, line=0.5, at=-3.8, font=2, cex=1.1)

#####################################################################################

dev.off()

#####################################################################################
