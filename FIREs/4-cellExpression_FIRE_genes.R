######## [[1]] Promoter
### [1] Processing the data
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
setwd("/proj/hyejunglab/NeuN/FIRE")
NeuNposgene = read.table("neuron_diffFIRE_promoter_gene.txt", header=F,sep="\t")[,1]  ##no extend
NeuNneggene = read.table("glia_diffFIRE_promoter_gene.txt", header=F,sep="\t")[,1]  ##no extend

targetgene = list(NeuNposgene, NeuNneggene)
diseasename = c("NeuNpos", "NeuNneg")

cellexp<-read.table("Capstone4_singlecellexpression_adultbrain.txt", header=T)
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

p.brainExp1 = ggplot(dat,aes(x=celltype, y=Expr, fill=Group)) + geom_hline(yintercept=0, colour="red") + #scale_fill_manual(values=c("#99CA3C","#75B3E1")) + 
  ylab("Normalized expression") +  geom_boxplot(outlier.size = NA) + theme(axis.text.x=element_text(angle = 90, hjust=1)) +
  ggtitle("Expression in each cell type: Promoter") 

######## [[2]] Enhancer-promoter interactions
### [1] Overlap FIREs with enhancer-promoter interactions based on Hi-C
NeuNposint = read.table("TSS_interacting_region_NeuNpos_FDR.bed")
NeuNnegint = read.table("TSS_interacting_region_NeuNneg_FDR.bed")

NeuNposint[NeuNposint[,1]=="chr23", 1] <- "chrX"
NeuNnegint[NeuNnegint[,1]=="chr23", 1] <- "chrX"

NeuNposregion = GRanges(seqnames=NeuNposint[,1], IRanges(NeuNposint[,4], NeuNposint[,4]+10000), chrint=NeuNposint[,2])
NeuNnegregion = GRanges(seqnames=NeuNnegint[,1], IRanges(NeuNnegint[,4], NeuNnegint[,4]+10000), chrint=NeuNnegint[,2])

olap = findOverlaps(NeuNposregion, promoteranges)
NeuNpospromoter = NeuNposregion[queryHits(olap)]
mcols(NeuNpospromoter) = cbind(mcols(NeuNposregion[queryHits(olap)]),mcols(promoteranges[subjectHits(olap)]))

olap = findOverlaps(NeuNnegregion, promoteranges)
NeuNnegpromoter = NeuNnegregion[queryHits(olap)]
mcols(NeuNnegpromoter) = cbind(mcols(NeuNnegregion[queryHits(olap)]),mcols(promoteranges[subjectHits(olap)]))

NeuNposenhancer = GRanges(seqnames=seqnames(NeuNpospromoter), IRanges(NeuNpospromoter$chrint, NeuNpospromoter$chrint+10000), tss=start(NeuNpospromoter), gene=NeuNpospromoter$gene)
NeuNnegenhancer = GRanges(seqnames=seqnames(NeuNnegpromoter), IRanges(NeuNnegpromoter$chrint, NeuNnegpromoter$chrint+10000), tss=start(NeuNnegpromoter), gene=NeuNnegpromoter$gene)

olap = findOverlaps(NeuNposenhancer, firanges.NeuNpos)
NeuNposfire = NeuNposenhancer[queryHits(olap)]

olap = findOverlaps(NeuNnegenhancer, firanges.NeuNneg)
NeuNnegfire = NeuNnegenhancer[queryHits(olap)]

### [2] Plot expression of genes mapped to FIREs
NeuNposensg = unique(NeuNposfire$gene)
NeuNnegensg = unique(NeuNnegfire$gene)

NeuNposgene = geneAnno1[geneAnno1$ensembl_gene_id %in% NeuNposensg, "hgnc_symbol"]
NeuNneggene = geneAnno1[geneAnno1$ensembl_gene_id %in% NeuNnegensg, "hgnc_symbol"]
NeuNposgene = NeuNposgene[NeuNposgene!=""]
NeuNneggene = NeuNneggene[NeuNneggene!=""]

##########[4] write the gene symbols to file
write.table(file="NeuNposgene_FIRE_loop.txt",as.data.frame(NeuNposgene),sep='\t',quote=F,row.names=F)
write.table(file="NeuNneggene_FIRE_loop.txt",as.data.frame(NeuNneggene),sep='\t',quote=F,row.names=F)


targetgene = list(NeuNposgene, NeuNneggene)
diseasename = c("NeuNpos", "NeuNneg")

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

p.brainExp2 = ggplot(dat,aes(x=celltype, y=Expr, fill=Group)) + geom_hline(yintercept=0, colour="red") + #scale_fill_manual(values=c("#99CA3C","#75B3E1")) + 
  ylab("Normalized expression") +  geom_boxplot(outlier.size = NA) + theme(axis.text.x=element_text(angle = 90, hjust=1)) +
  ggtitle("Expression in each cell type: Hi-C") 

pdf("FIRE_matched_gene_cell-type_expression.pdf", width=5, height=5)
p.brainExp1
p.brainExp2
dev.off()

