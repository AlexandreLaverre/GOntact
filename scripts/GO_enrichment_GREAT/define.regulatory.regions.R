#######################################################################

pathEnsembl="../../data/ensembl_annotations/"
pathGO="../../data/GeneOntology/"
pathResults="../../results/GO_enrichment_GREAT/"

ensrelease=94

chromo=c(as.character(1:22), "X", "Y")

#######################################################################

pars=list()
pars[["classical_upstream5kb_downstream1kb_extend100kb"]]=list("basal5"=5e3, "basal3"=1e3, "max.extend"=1e5, "min.dist"=0)
pars[["classical_upstream5kb_downstream1kb_extend1Mb"]]=list("basal5"=5e3, "basal3"=1e3, "max.extend"=1e6, "min.dist"=0)
pars[["classical_upstream0kb_downstream0kb_extend1Mb"]]=list("basal5"=0, "basal3"=0, "max.extend"=1e6, "min.dist"=0)
pars[["classical_mindist25kb_extend1Mb"]]=list("basal5"=0, "basal3"=0, "max.extend"=1e6, "min.dist"=25000)
pars[["classical_mindist25kb_extend2Mb"]]=list("basal5"=0, "basal3"=0, "max.extend"=2e6, "min.dist"=25000)

#######################################################################

general.GO=c("GO:0005575", "GO:0110165", "GO:0003674", "GO:0008150")

## cellular component
## cellular anatomical entity
## molecular function
## biological process

#######################################################################

for(space in c("biological_process", "molecular_function", "cellular_component")){

  for(method in names(pars)){
    basal5=pars[[method]][["basal5"]]
    basal3=pars[[method]][["basal3"]]
    min.dist=pars[[method]][["min.dist"]]
    max.extend=pars[[method]][["max.extend"]]
    outputdir=method
    
    for(sp in c("human", "mouse")){
      ## check if already done

      pathOut=paste(pathResults, sp, "/",outputdir,"/regulatory_regions_Ensembl",ensrelease,"_",space,".txt",sep="")
      
      if(file.exists(pathOut)){
        print(paste("already done", method, sp))
      } else{
        dirOut=paste(pathResults, sp, "/",outputdir,sep="")
        
        if(!dir.exists(dirOut)){
          system(paste("mkdir ",pathResults, sp, "/",outputdir,sep=""))
        }
        
        ## chr sizes
        
        chr.sizes=read.table(paste(pathEnsembl, sp, "/chr_sizes_Ensembl",ensrelease,".txt",sep=""),h=F, stringsAsFactors=F)
        colnames(chr.sizes)=c("chr", "size")
        
        ## gene info
        gene.info=read.table(paste(pathEnsembl, sp, "/GeneInfo_Ensembl",ensrelease,".txt",sep=""), h=T, stringsAsFactors=F, sep="\t", quote="\"")
        
        print(paste(nrow(gene.info), "genes"))
        
        ## protein-coding genes
        gene.info=gene.info[which(gene.info$biotype=="protein_coding"),]
        
        print(paste(nrow(gene.info), "pc genes"))
        
        ## assembled chromosomes
        gene.info=gene.info[which(gene.info$name%in%chromo),]
        
        print(paste(nrow(gene.info), "genes on assembled chromo"))
        
        ## gene names
        gene.names=read.table(paste(pathEnsembl, sp, "/GeneNames_Ensembl",ensrelease,".txt",sep=""), h=T, stringsAsFactors=F, sep="\t", quote="\"")
        gene.names=gene.names[which(gene.names$Gene.stable.ID%in%gene.info$stable_id),]
        rownames(gene.names)=gene.names$Gene.stable.ID
        
        gene.info$Symbol=gene.names[gene.info$stable_id,"Gene.name"]
        
        ## gene ontology
        ## after removing GO annotations
        go=read.table(paste(pathGO, sp, ".simplified.gene.annotation.",space,".txt",sep=""),h=T, sep="\t", stringsAsFactors=F, quote="\"")
        print(paste(length(unique(go$GeneName)), "genes in GO annotation"))
        
        go=go[which(!go$GOID%in%general.GO),]
        
        print(paste(length(unique(go$GeneName)), "genes in GO annotation after removing non-meaningful categories"))
        
        ## select only genes that have at least one GO category
        gene.info=gene.info[which(gene.info$Symbol%in%go$GeneName),]
        
        print(paste(nrow(gene.info), "genes with GO annotations"))
        
        ## canonical transcript TSS
        
        tss=read.table(paste(pathEnsembl, sp, "/canonical_transcripts_Ensembl",ensrelease,"_protein_coding.txt",sep=""), h=T, stringsAsFactors=F, sep="\t", quote="\"")
        tss=tss[which(tss$gene_id%in%gene.info$stable_id),]
        
        print(paste(nrow(tss), "genes with canonical TSS for pc transcripts"))
        
        ## regulatory regions
        
        reg.regions=tss[,c("gene_id", "chr", "strand", "TSS")]
        reg.regions$gene_name=gene.names[reg.regions$gene_id, "Gene.name"]
        
        ## reorder columns
        reg.regions=reg.regions[,c("gene_id","gene_name", "chr", "strand", "TSS")]
        rownames(reg.regions)=reg.regions$gene_id
        
        reg.regions$region_start=rep(NA, nrow(reg.regions))
        reg.regions$region_end=rep(NA, nrow(reg.regions))
        
        reg.regions$exclude_start=rep(NA, nrow(reg.regions))
        reg.regions$exclude_end=rep(NA, nrow(reg.regions))
        
        for(chr in chr.sizes$chr){
          print(paste("chr", chr))
          
          this.chrsize=chr.sizes$size[which(chr.sizes$chr==chr)]
          
          this.tss=tss[which(tss$chr==chr),]
          this.tss=this.tss[order(this.tss$TSS),] ## we sort TSS
          
          for(i in 1:nrow(this.tss)){
            this.gene=this.tss$gene_id[i]
            
            this.strand=this.tss$strand[i]
            
            this.start=NA
            this.end=NA
            
            if(this.strand==1){
              this.start=this.tss$TSS[i]-basal5
              this.end=this.tss$TSS[i]+basal3-1
            }
            
            if(this.strand==-1){
              this.start=this.tss$TSS[i]-basal3+1
              this.end=this.tss$TSS[i]+basal5
            }
            
            ## check chromosome
            
            if(this.start<1){
              this.start=1
            }
            
            if(this.end>this.chrsize){
              this.end=this.chrsize
            }
            
            ## print(paste(this.gene, "basal", this.start, this.end))
            
            ## extend left
            
            if(i>1){
              previous.strand=this.tss$strand[i-1]
              previous.basal=NA
              
              if(previous.strand==1){
                previous.basal=this.tss$TSS[i-1]+basal3-1
              }
              
              if(previous.strand==-1){
                previous.basal=this.tss$TSS[i-1]+basal5
              }
              
              this.minstart=max(c(1, previous.basal+1, this.tss$TSS[i]-max.extend))
              
              ##  print(paste("left basal", previous.basal))
              
              if(this.minstart<this.start){
                this.start=this.minstart
              }
            } else{
              ## no gene before, we extend till the start of the chromosome
              this.minstart=max(1, this.tss$TSS[i]-max.extend)
              
              if(this.minstart<this.start){
                this.start=this.minstart
              }
            }
            
            ## extend right
            
            if(i<nrow(this.tss)){
              next.strand=this.tss$strand[i+1]
              next.basal=NA
              
              if(next.strand==1){
                next.basal=this.tss$TSS[i+1]-basal5
              }
              
              if(next.strand==-1){
                next.basal=this.tss$TSS[i+1]-basal3+1
              }
              
              ## print(paste("right basal", next.basal))
              
              this.maxend=min(this.chrsize, next.basal-1, this.tss$TSS[i]+max.extend-1)
              
              if(this.maxend>this.end){
                this.end=this.maxend
              }
            } else{
              ## no other gene, end of chr
              this.maxend=min(this.chrsize, this.tss$TSS[i]+max.extend-1)
              
              if(this.maxend>this.end){
                this.end=this.maxend
              }
            }
            
            ##  print(paste("final", this.start, this.end))
            
            ## we exclude a region around in the TSS
            
            if(min.dist>0){
              reg.regions[this.gene,  "exclude_start"]=this.tss$TSS[i]-min.dist
              reg.regions[this.gene,  "exclude_end"]=this.tss$TSS[i]+min.dist
            }
            
            ## fill in info
            
            reg.regions[this.gene,  "region_start"]=this.start
            reg.regions[this.gene,  "region_end"]=this.end
            
          }
        }
        
        write.table(reg.regions, file=pathOut, row.names=F, col.names=T, sep="\t", quote=F)
        
      }
    }
  }
}

#######################################################################
