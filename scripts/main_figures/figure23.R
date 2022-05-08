#####################################################################################

pathFigures="../../data_for_publication/figures/"

options(stringsAsFactors=F)

#####################################################################################

if(load){
  load(paste(pathFigures, "RData/data.enrichment.results.RData",sep=""))

  species=c("human", "human")
  samples=c("Vista_heart_vs_ENCODE", "Vista_midbrain_vs_ENCODE")
  domain="biological_process"

  minFDR=1e-10
  maxFDR=0.05
  minEnrichment=1.5

  load=FALSE
}

#####################################################################################

if(prepare){

  results=list()

  methods=c("GREAT_upstream5kb_downstream1kb_extend1Mb", "GREAT_fixed_size_upstream1Mb_downstream1Mb", "contacts_mindist0kb_maxdist1Mb",  "shared_contacts_mindist0kb_maxdist1Mb", "hybrid_mindist25kb_maxdist1Mb")
  shortnames=c("GREAT (basal + extension)", "fixed 1Mb window", "GOntact (all PCHi-C data)", "GOntact (common contacts)", "GOntact (hybrid)")

  for(i in 1:2){
    sp=species[i]
    sample=samples[i]

    results[[paste(sp, sample)]]=list()

    for(j in 1:length(methods)){
      method=methods[j]
      name=shortnames[j]

      ## enrichment results

      enrichment=enrichment.results[[sp]][[sample]][[domain]][[method]]
      rownames(enrichment)=enrichment$GOID
      enrichment$Enrichment=enrichment$Observed/enrichment$Expected
      enrichment$FDR[which(enrichment$FDR<minFDR)]=minFDR  ## replace very small FDRs
      signif=enrichment$GOID[which(enrichment$FDR<maxFDR  & enrichment$Enrichment>=minEnrichment)]

      ## smaller names

      enrichment$GOName[which(enrichment$GOName=="positive regulation of epithelial cell proliferation")]="positive reg. epithelial cell proliferation"
      enrichment$GOName[which(enrichment$GOName=="regulation of cellular macromolecule biosynthetic process")]="macromolecule biosynthetic process"
      enrichment$GOName[which(enrichment$GOName=="positive regulation of transcription by RNA polymerase II")]="positive reg. of transcription by RNA pol II"
      enrichment$GOName[which(enrichment$GOName=="regulation of nucleobase-containing compound metabolic process")]="reg. nucleobase-containing compound metabolic process"

      ## save results

      results[[paste(sp, sample)]][[name]]=list("enrichment"=enrichment, "signif"=signif)

    }
  }

  prepare=FALSE
}

#####################################################################################

## 1 column width 85 mm = 3.34 in
## 1.5 column width 114 mm = 4.49 in
## 2 columns width 174 mm = 6.85 in
## max height: 11 in

pdf(file=paste(pathFigures, "Figure2.pdf",sep=""), width=6.85, height=6)

m=matrix(rep(NA, 12*16), nrow=12)

for(i in 1:6){
  m[i,]=c(rep(1,4), rep(2,4), rep(5,4), rep(6,4))
}

for(i in 7:12){
  m[i,]=c(rep(3,4), rep(4,4), rep(7,4), rep(8,4))
}


layout(m)

par(oma=c(0,1,1.5,0))

#####################################################################################

labels=toupper(letters[1:8])

sp="human"

i=1

label.pos=c(0.25,0.75)
names(label.pos)=c("heart", "midbrain")

for(sample in c("Vista_heart_vs_ENCODE", "Vista_midbrain_vs_ENCODE")){
  tissue=unlist(strsplit(sample, split="_"))[2]

  for(method in shortnames[1:2]){

    enrichment=results[[paste(sp, sample)]][[method]][["enrichment"]]
    signif=results[[paste(sp, sample)]][[method]][["signif"]]

    top10=signif[10:1]

    ## barplot, observed vs expected
    observed=enrichment[top10,"Observed"]
    expected=enrichment[top10,"Expected"]

    all.values=-100*as.numeric(rbind(expected, observed))
    xax=pretty(all.values)

    xlim=range(c(xax, all.values))

    par(mar=c(3.5, 1.5, 1, 0.35))
    b=barplot(all.values, horiz=T, axes=F, col=c("black","indianred"), border=NA, xlim=xlim, space=c(2.25,0.35))

    axis(side=1, mgp=c(3,0.5,0), at=xax, labels=-xax, cex.axis=0.95)

    mtext(labels[i], side=3, at=xlim[1]-diff(xlim)/10, line=0.15, font=2, cex=1.05)
    mtext("% associated enhancers", side=1, line=1.75, cex=0.75)

    if(i==8){
      legend("bottomleft", fill=c("indianred", "black"), border=c("indianred", "black"), legend=c("observed", "expected"), inset=c(-0.25, 0.02), cex=0.95, xpd=NA)
    }

    if(sample=="Vista_heart_vs_ENCODE"){
      mtext(method, side=2, cex=0.8, font=1, line=0.5)
    }
    ## barplot, log10 FDR

    log10fdr=-log10(enrichment[top10,"FDR"])

    par(mar=c(3.5, 0.35, 1, 1.5))
    b=barplot(log10fdr, horiz=T, axes=F, col="steelblue", border=NA, space=1, xlim=c(0,10))

    xax=seq(from=0, to=10, by=2)
    xax.labels=as.character(xax)
    xax.labels[length(xax.labels)]=">10"
    axis(side=1, mgp=c(3,0.5,0), cex.axis=0.95, at=xax, labels=xax.labels)

    mtext("-log10 FDR", side=1, line=1.75, cex=0.75)

    text(enrichment[top10, "GOName"], y=b+1, x=0, adj=c(0.5,0.5), cex=0.95, xpd=NA)

    i=i+1

  }

  mtext(paste("Vista",tissue), side=3, line=0.25, outer=T, cex=0.85, at=label.pos[tissue], font=2)

}

dev.off()

#####################################################################################
#####################################################################################

pdf(file=paste(pathFigures, "Figure3.pdf",sep=""), width=6.85, height=8.2)

m=matrix(rep(NA, 18*16), nrow=18)

for(i in 1:6){
  m[i,]=c(rep(1,4), rep(2,4), rep(7,4), rep(8,4))
}

for(i in 7:12){
  m[i,]=c(rep(3,4), rep(4,4), rep(9,4), rep(10,4))
}

for(i in 13:18){
  m[i,]=c(rep(5,4), rep(6,4), rep(11,4), rep(12,4))
}

layout(m)

par(oma=c(0,1,1.5,0))

#####################################################################################

labels=toupper(letters[1:8])

sp="human"

i=1

label.pos=c(0.25,0.75)
names(label.pos)=c("heart", "midbrain")

for(sample in c("Vista_heart_vs_ENCODE", "Vista_midbrain_vs_ENCODE")){
  tissue=unlist(strsplit(sample, split="_"))[2]

  for(method in shortnames[3:5]){

    enrichment=results[[paste(sp, sample)]][[method]][["enrichment"]]
    signif=results[[paste(sp, sample)]][[method]][["signif"]]

    top10=signif[10:1]

    ## barplot, observed vs expected
    observed=enrichment[top10,"Observed"]
    expected=enrichment[top10,"Expected"]

    all.values=-100*as.numeric(rbind(expected, observed))
    xax=pretty(all.values)

    xlim=range(c(xax, all.values))

    par(mar=c(3.5, 1.5, 1, 0.35))
    b=barplot(all.values, horiz=T, axes=F, col=c("black","indianred"), border=NA, xlim=xlim, space=c(2.25,0.35))

    axis(side=1, mgp=c(3,0.5,0), at=xax, labels=-xax, cex.axis=0.95)

    mtext(labels[i], side=3, at=xlim[1]-diff(xlim)/10, line=0.15, font=2, cex=1.05)
    mtext("% associated enhancers", side=1, line=1.75, cex=0.75)

    if(i==8){
      legend("bottomleft", fill=c("indianred", "black"), border=c("indianred", "black"), legend=c("observed", "expected"), inset=c(-0.25, 0.02), cex=0.95, xpd=NA)
    }

    if(sample=="Vista_heart_vs_ENCODE"){
      mtext(method, side=2, cex=0.8, font=1, line=0.5)
    }
    ## barplot, log10 FDR

    log10fdr=-log10(enrichment[top10,"FDR"])

    par(mar=c(3.5, 0.35, 1, 1.5))
    b=barplot(log10fdr, horiz=T, axes=F, col="steelblue", border=NA, space=1, xlim=c(0,10))

    xax=seq(from=0, to=10, by=2)
    xax.labels=as.character(xax)
    xax.labels[length(xax.labels)]=">10"
    axis(side=1, mgp=c(3,0.5,0), cex.axis=0.95, at=xax, labels=xax.labels)

    mtext("-log10 FDR", side=1, line=1.75, cex=0.75)

    text(enrichment[top10, "GOName"], y=b+1, x=0, adj=c(0.5,0.5), cex=0.95, xpd=NA)

    i=i+1

  }

  mtext(paste("Vista",tissue), side=3, line=0.25, outer=T, cex=0.85, at=label.pos[tissue], font=2)

}

dev.off()

#####################################################################################
#####################################################################################
