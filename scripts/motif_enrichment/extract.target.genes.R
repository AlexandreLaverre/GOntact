################################################################################

pathTargetGenes="../../data/putative_target_genes/"
pathResults="../../results/"

sp="mouse"
sample="Vista_heart_vs_ENCODE"

options(stringsAsFactors=F)

################################################################################

## great
great=read.table(paste(pathResults, sp, "/",sample, "/biological_process/GREAT_upstream5kb_downstream1kb_extend1Mb/GOntact_element_gene_association_background.txt", sep=""), h=T, stringsAsFactors=F)

## contacts
gontact=read.table(paste(pathResults, sp, "/", sample, "/biological_process/contacts_mindist0kb_maxdist1Mb/GOntact_element_gene_association_background.txt",sep=""), h=T, stringsAsFactors=F)

################################################################################

hnf4a=read.table(paste(pathTargetGenes, "Marable2020_HNF4a_knockout.txt",sep=""), h=T, stringsAsFactors=F, sep="\t", quote="\"", comment.char="")

down=hnf4a$gene[which(hnf4a$log2.fold_change.<0 & hnf4a$q_value<0.05)]

################################################################################

elements.gontact=unique(gontact$ElementID[which(gontact$GeneSymbol%in%down)])
elements.great=unique(great$ElementID[which(great$GeneSymbol%in%down)])

only.gontact=setdiff(elements.gontact, elements.great)
only.great=setdiff(elements.great, elements.gontact)
shared=intersect(elements.great, elements.gontact)

################################################################################

outdir="HNF4a"

system(paste("mkdir -p ",pathResults, "motif_enrichment/",sp,"/",outdir,"/GREAT_only/", sep=""))
system(paste("mkdir -p ",pathResults, "motif_enrichment/",sp,"/",outdir,"/GOntact_only/", sep=""))
system(paste("mkdir -p ",pathResults, "motif_enrichment/",sp,"/",outdir,"/shared/", sep=""))


## write output for only great
chr.great=unlist(lapply(only.great, function(x) unlist(strsplit(x, split=":"))[1]))
coords.great=unlist(lapply(only.great, function(x) unlist(strsplit(x, split=":"))[2]))
start.great=unlist(lapply(coords.great, function(x) unlist(strsplit(x, split="-"))[1]))
end.great=unlist(lapply(coords.great, function(x) unlist(strsplit(x, split="-"))[2]))

res.great=data.frame(chr.great, start.great, end.great, only.great, rep(".", length(only.great)),rep("+", length(only.great)))

write.table(res.great, paste(pathResults, "motif_enrichment/",sp,"/",outdir,"/GREAT_only/associated_enhancers.bed",sep=""), row.names=F, col.names=F, sep="\t", quote=F)

## write output for only gontact

chr.gontact=unlist(lapply(only.gontact, function(x) unlist(strsplit(x, split=":"))[1]))
coords.gontact=unlist(lapply(only.gontact, function(x) unlist(strsplit(x, split=":"))[2]))
start.gontact=unlist(lapply(coords.gontact, function(x) unlist(strsplit(x, split="-"))[1]))
end.gontact=unlist(lapply(coords.gontact, function(x) unlist(strsplit(x, split="-"))[2]))

res.gontact=data.frame(chr.gontact, start.gontact, end.gontact, only.gontact, rep(".", length(only.gontact)),rep("+", length(only.gontact)))

write.table(res.gontact, paste(pathResults, "motif_enrichment/",sp,"/",outdir,"/GOntact_only/associated_enhancers.bed",sep=""), row.names=F, col.names=F, sep="\t", quote=F)

## write output for shared

chr.shared=unlist(lapply(shared, function(x) unlist(strsplit(x, split=":"))[1]))
coords.shared=unlist(lapply(shared, function(x) unlist(strsplit(x, split=":"))[2]))
start.shared=unlist(lapply(coords.shared, function(x) unlist(strsplit(x, split="-"))[1]))
end.shared=unlist(lapply(coords.shared, function(x) unlist(strsplit(x, split="-"))[2]))

res.shared=data.frame(chr.shared, start.shared, end.shared, shared, rep(".", length(shared)),rep("+", length(shared)))

write.table(res.shared, paste(pathResults, "motif_enrichment/",sp,"/",outdir,"/shared/associated_enhancers.bed",sep=""), row.names=F, col.names=F, sep="\t", quote=F)

################################################################################
