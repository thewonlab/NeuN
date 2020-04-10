options(stringsAsFactors=F)
library(stats4)
library(BiocGenerics)
library(parallel)
library(S4Vectors)
library(IRanges)
library(GenomeInfoDb)
library(GenomicRanges)
library(ggplot2)
library(reshape2)
library(dplyr)
#setwd("/proj/hyejunglab/NeuN/FIRE/")
load("geneAnno_allgenes.rda") # the saved file is geneAnno1
#ab<-load("geneAnno_allgenes.rda")
### [1] Processing the data
promoter = read.table("gencode19_promoter.bed")
cellexp = read.table("Capstone4_singlecellexpression_adultbrain.txt", header=T)
fire.NeuNneg = read.table("super_FIRE_call_NeuNneg.txt", header=T,sep="\t")
fire.NeuNpos = read.table("super_FIRE_call_NeuNpos.txt", header=T,sep="\t")

firanges.NeuNpos = GRanges(seqnames=fire.NeuNpos$chr, IRanges(fire.NeuNpos$start, fire.NeuNpos$end))
firanges.NeuNneg = GRanges(seqnames=fire.NeuNneg$chr, IRanges(fire.NeuNneg$start, fire.NeuNneg$end))

promoteranges = GRanges(seqnames=promoter[,1], IRanges(promoter[,2],promoter[,3]), gene=promoter[,6])

### [2] Overlap
olap = findOverlaps(firanges.NeuNpos, promoteranges)
firepro.NeuNpos = firanges.NeuNpos[queryHits(olap)]
mcols(firepro.NeuNpos) = cbind(mcols(firepro.NeuNpos), mcols(promoteranges[subjectHits(olap)]))

olap = findOverlaps(firanges.NeuNneg, promoteranges)
firepro.NeuNneg = firanges.NeuNneg[queryHits(olap)]
mcols(firepro.NeuNneg) = cbind(mcols(firepro.NeuNneg), mcols(promoteranges[subjectHits(olap)]))

### [3] Plot expression 
NeuNposensg = unique(firepro.NeuNpos$gene)
NeuNnegensg = unique(firepro.NeuNneg$gene)

NeuNposgene = geneAnno1[geneAnno1$ensembl_gene_id %in% NeuNposensg, "hgnc_symbol"]
NeuNneggene = geneAnno1[geneAnno1$ensembl_gene_id %in% NeuNnegensg, "hgnc_symbol"]

NeuNposgene = NeuNposgene[NeuNposgene!=""]
NeuNneggene = NeuNneggene[NeuNneggene!=""]

##########[4] write the gene symbols to file
write.table(file="NeuNposgene_superFIRE_H3K27ac.txt",as.data.frame(NeuNposgene),sep='\t',quote=F,row.names=F)
write.table(file="NeuNneggene_superFIRE_H3K27ac.txt",as.data.frame(NeuNneggene),sep='\t',quote=F,row.names=F)



targetgene = list(NeuNposgene, NeuNneggene)
diseasename = c("NeuNpos", "NeuNneg")

rownames(cellexp) = cellexp[,1]
cellexp = cellexp[,-1]

datExpr = scale(cellexp,center=T, scale=F) 

exprdat = vector(mode="list", length=length(targetgene))

for(i in 1:length(targetgene)){
  exprdat[[i]] = apply(datExpr[match(targetgene[[i]], rownames(datExpr)),],2,mean,na.rm=T)
}

dat = c()

for(i in 1:length(targetgene)){
  datframe = data.frame(Group=diseasename[i], cell=names(exprdat[[i]]), Expr=exprdat[[i]])
  dat = rbind(dat, datframe)
}

dat$cell = gsub("Fetal.quiescent","Fetal-quiescent",dat$cell)
dat$cell = gsub("Fetal.replicating","Fetal-replicating",dat$cell)

dat$celltype = unlist(lapply(strsplit(dat$cell, split="[.]"),'[[',1))

dat = dat[!(dat$celltype %in% c("OPC","Fetal-quiescent","Fetal-replicating",
                              "Ex1","Ex2","Ex3","Ex4","Ex5","Ex6","Ex7","Ex8",
                              "In1","In2","In3","In4","In5","In6","In7","In8")), ]

dat$celltype = factor(dat$celltype, levels=c("Neurons","Astrocytes","Microglia","Endothelial","Oligodendrocytes"))

pdf("superFIRE_matched_gene_cell-type_expression.pdf")

p.brainExp1 = ggplot(dat,aes(x=celltype, y=Expr, fill=Group)) + geom_hline(yintercept=0, colour="red") + #scale_fill_manual(values=c("#99CA3C","#75B3E1")) + 
  ylab("Normalized expression") +  geom_boxplot(outlier.size = NA) + theme(axis.text.x=element_text(angle = 90, hjust=1)) +
  ggtitle("superFIRE-associated gene\n 
          expression in each cell type: Promoter") 

p.brainExp1
dev.off()