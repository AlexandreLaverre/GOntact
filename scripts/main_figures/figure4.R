####################################################################################

pathFigures="../../data_for_publication/figures/"

options(stringsAsFactors=F)

#####################################################################################

if(load){
  load(paste(pathFigures, "RData/data.enrichment.results.RData",sep=""))

  sp="human"
  samples=c("Vista_heart_vs_ENCODE", "Vista_midbrain_vs_ENCODE")
  domain="biological_process"

  minFDR=1e-10
  maxFDR=0.05
  minEnrichment=1.5

  load=FALSE
}

#####################################################################################

if(prepare){

  signif.list=list()
  all.signif=list()
  presence=list()

  methods=c("GREAT_upstream5kb_downstream1kb_extend1Mb", "GREAT_fixed_size_upstream1Mb_downstream1Mb", "contacts_mindist0kb_maxdist1Mb",  "shared_contacts_mindist0kb_maxdist1Mb", "hybrid_mindist25kb_maxdist1Mb")
  shortnames=c("GREAT\n(basal + extension)", "fixed 2Mb window", "GOntact\n(all PCHi-C data)", "GOntact\n(common contacts)", "GOntact\n(hybrid)")
  names(shortnames)=methods

  col.methods=c("red", "black", "darkorange", "slateblue", "steelblue")
  names(col.methods)=methods


  prop.expected=list()

  for(sample in samples){
    signif.list[[sample]]=list()
    prop.expected[[sample]]=list()

    for(method in methods){
      test=enrichment.results[[sp]][[sample]][[domain]][[method]]
      test$Enrichment=test$Observed/test$Expected
      rownames(test)=test$GOID
      signif=test$GOID[which(test$FDR < maxFDR & test$Enrichment >= minEnrichment)]

      signif.list[[sample]][[method]]=signif

      prop.expected[[sample]][[method]]=test[signif, "Expected"]
    }

    this.all.signif=unique(unlist(signif.list[[sample]]))

    this.presence=matrix(rep(NA, length(this.all.signif)*length(methods)), nrow=length(this.all.signif))
    colnames(this.presence)=methods

    for(method in methods){
      this.presence[,method]=as.numeric(this.all.signif%in%signif.list[[sample]][[method]])
    }

    presence[[sample]]=this.presence

    all.signif[[sample]]=this.all.signif
  }

  ## combinations

  possibilities=list()
  for(method in methods){
    possibilities[[method]]=c(0,1)
  }

  combinations=list()

  for(i1 in possibilities[[methods[1]]]){
    for(i2 in possibilities[[methods[2]]]){
      for(i3 in possibilities[[methods[3]]]){
        for(i4 in possibilities[[methods[4]]]){
          for(i5 in possibilities[[methods[5]]]){
            n=length(combinations)
            combinations[[n+1]]=c(i1,i2,i3,i4,i5)
          }
        }
      }
    }
  }

  combinations=as.data.frame(combinations)
  rownames(combinations)=methods
  colnames(combinations)=1:ncol(combinations)

  ## values

  values=list()

  for(sample in samples){
    this.values=rep(NA, ncol(combinations))
    this.presence=presence[[sample]]

    for(i in 1:ncol(combinations)){
      this.comb=combinations[,i]
      this.values[i]= length(which(unlist(lapply(1:nrow(this.presence), function(x) all(as.numeric(this.presence[x,])==as.numeric(this.comb))))))
    }

    values[[sample]]=this.values
  }


  prepare=FALSE
}

#####################################################################################

## 1 column width 85 mm = 3.34 in
## 1.5 column width 114 mm = 4.49 in
## 2 columns width 174 mm = 6.85 in
## max height: 11 in


pdf(file=paste(pathFigures, "Figure4.pdf",sep=""), width=6.85, height=9)

m=matrix(rep(NA, 27*10), nrow=27)

for(i in 1:6){
  m[i,]=c(rep(1, 10))
}

for(i in 7:10){
  m[i,]=c(rep(2, 10))
}


for(i in 11:16){
  m[i,]=c(rep(4, 10))
}

for(i in 17:20){
  m[i,]=c(rep(5, 10))
}


for(i in 21:27){
  m[i,]=c(rep(3, 5), rep(6,5))
}


layout(m)

#####################################################################################

## upset plot - barplot for the values

sample.names=c("Vista heart", "Vista midbrain")
names(sample.names)=samples


labels=list()
labels[[samples[1]]]=c("A", "C")
labels[[samples[2]]]=c("B", "D")

for(sample in samples){

  this.values=values[[sample]]
  ylim=c(0, max(this.values))

  xpos=1:length(this.values)
  xlim=c(0.5, length(this.values)+0.5)

  par(mar=c(0.5,7.5, 2.5,2.5))
  plot(1, type="n", xlab="", ylab="", xlim=xlim, ylim=ylim, yaxs="i", xaxs="i", axes=F)

  rect(xpos-0.25, 0, xpos+0.25, this.values, col="black", border=NA)

  axis(side=2, mgp=c(3, 0.65, 0), cex.axis=0.95, las=2)
  axis(side=1, at=xpos, labels=rep("", length(xpos)), col="gray60")

  mtext("nb. significant GO categories", side=2, line=2.75, cex=0.75)

  mtext(sample.names[sample], side=3, cex=0.85, line=0.5)

  mtext(labels[[sample]][1], at=-2, font=2, cex=1.1, line=1)

 #####################################################################################

  ## labels for the intersections

  nb.methods=length(methods)
  ylim=c(0.5,nb.methods+0.5)
  xpos=1:length(this.values)
  xlim=c(0.5, max(xpos)+0.5)
  ypos=nb.methods:1
  names(ypos)=methods

  par(mar=c(0.1,7.75,0.5,2.5))
  plot(1, type="n", xlab="", ylab="", xlim=xlim, ylim=ylim, xaxs="i", yaxs="i", axes=F)
  smallx=diff(xlim)/40

  mtext(shortnames, side=2,  at=ypos, las=2, cex=0.65, line=0.5, col=col.methods)

  for(method in methods){
    segments(xlim[1], ypos[method], xlim[2], ypos[method], lty=3)

    drawx=xpos[which(combinations[method,]==1)]
    points(drawx, rep(ypos[method], length(drawx)), cex=1.1, pch=20, col=col.methods[method])

    mtext(length(signif.list[[sample]][[method]]), side=4, at=ypos[method], cex=0.6, las=2)
  }

  ##expected proportion for all methods

  xpos.expected=1:5
  names(xpos.expected)=methods

  if(sample=="Vista_heart_vs_ENCODE"){
    ylim.expected=c(0, 25)
  }

  if(sample=="Vista_midbrain_vs_ENCODE"){
    ylim.expected=c(0,10)
  }

  xlim.expected=c(0.25, 5.75)


  par(mar=c(1.1, 5.1, 2.1, 1.1))
  plot(1, type="n", xlab="", ylab="", axes=F, xlim=xlim.expected, ylim=ylim.expected, xaxs="i")

  for(method in methods){
    boxplot(100*prop.expected[[sample]][[method]], add=T, at=xpos.expected[method], col="white", border=col.methods[method], axes=F, boxwex=0.85, outline=F, notch=F)
  }

  axis(side=2, mgp=c(3,0.5,0), cex.axis=0.95)

  axis(side=1, at=xpos.expected, labels=rep("", length(methods)))

  mtext("% associated ENCODE enhancers\nsignificant GO categories ", side=2, line=2.5, cex=0.75)

  mtext(labels[[sample]][2], at=-0.5, font=2, cex=1.1, line=0.5)
  mtext(sample.names[sample], side=3, cex=0.85, line=0.5)

}

#####################################################################################

dev.off()

#####################################################################################
