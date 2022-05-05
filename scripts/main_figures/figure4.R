#####################################################################################

pathFigures="../../data_for_publication/figures/"
pathJASPAR="../../data/JASPAR/"

options(stringsAsFactors=F)

library(imager)

#####################################################################################

if(load){

  load(paste(pathFigures, "RData/data.motif.enrichment.RData",sep=""))

  ## images

  images=list()
  images[["MA0909.3"]]=load.image(paste(pathJASPAR, "MA0909.3.png",sep=""))
  images[["MA0114.3"]]=load.image(paste(pathJASPAR, "MA0114.3.png",sep=""))
  images[["MA1537.1"]]=load.image(paste(pathJASPAR, "MA1537.1.png",sep=""))

  load=FALSE
}

#####################################################################################
#####################################################################################

## 1 column width 85 mm = 3.34 in
## 1.5 column width 114 mm = 4.49 in
## 2 columns width 174 mm = 6.85 in
## max height: 11 in

pdf(file=paste(pathFigures, "Figure4.pdf",sep=""), width=6.85, height=3.5)

m=matrix(rep(NA, 10*6), nrow=10)

for(i in 1:9){
  m[i,]=c(rep(1,2), rep(3,2), rep(5,2))
}

for(i in 10){
  m[i,]=c(rep(2,2), rep(4,2), rep(6,2))
}

layout(m)

#####################################################################################

labels=c("A", "B", "C")
names(labels)=names(images)

for(motif in names(images)){

  ## proportions with motif

  prop.with=100*nb.with[[motif]]/nb.tot[[motif]]

  ci.low=100*unlist(lapply(names(nb.with[[motif]]), function(x) prop.test(nb.with[[motif]][x], nb.tot[[motif]][x])$conf.int[1]))
  ci.high=100*unlist(lapply(names(nb.with[[motif]]), function(x) prop.test(nb.with[[motif]][x], nb.tot[[motif]][x])$conf.int[2]))

  ylim=c(0, max(ci.high)+1)

  par(mar=c(2.25, 4.1, 2.15, 0.5))
  b=barplot(prop.with, space=c(0.5, 1, 0.5, 0.5,  1), ylim=ylim, names=rep("",5), col=c("gray40", "black", "darkorange", "slateblue", "lightblue"), density=30, border=c("gray40", "black", "darkorange", "slateblue", "lightblue"))

  segments(b, ci.low, b,  ci.high, xpd=NA, lwd=1.5)
  tinyx=0.1
  segments(b-tinyx, ci.low, b+tinyx, ci.low)
  segments(b-tinyx, ci.high, b+tinyx, ci.high)

  ## logo

  img=images[[motif]]

  par(mar=c(0.5, 1, 0.15, 2.25))
  plot(img, axes=F, ylim=c(410, 40), xlim=c(0, 900))



}

#####################################################################################

dev.off()

#####################################################################################
