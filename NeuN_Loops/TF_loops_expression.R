options(stringsAsFactors=F)
library(ggplot2)
library(GenomicRanges)
library(plyr)
library(reshape)
setwd("/nas/longleaf/home/hubenxia/project/RNA-seq/neuron_vs_glia/known_gene/")
load("/proj/hyejunglab/NeuN/RNAseq/analysis/datExprMeta.rda")
FPKM<-read.table("neuron_glia_FPKM_dot_remove.txt",header=T,sep="\t",row.names=1)
colnames(FPKM)<-c("CMC_HBCC_028_neg","CMC_HBCC_025_pos","CMC_HBCC_028_pos","CMC_HBCC_027_pos",
                        "CMC_HBCC_025_neg","CMC_HBCC_027_neg","CMC_HBCC_016_pos","CMC_HBCC_016_neg")

FPKM_neuron<-as.data.frame((apply(FPKM[,c(2:4,7)],1,mean)))
FPKM_glia<-as.data.frame((apply(FPKM[,c(1,5:6,8)],1,mean)))
FPKM<-cbind(FPKM_neuron,FPKM_glia)
FPKM$gene<-rownames(FPKM)
rownames(FPKM)<-c()
colnames(FPKM)<-c('neuron','glia','geneid')
setwd("/nas/longleaf/home/hubenxia/project/ATAC-seq")
neuron<-read.table("neuron_TF_gene.txt",sep="\t",header=F)[,-3]
glia<-read.table("glia_TF_gene.txt",sep="\t",header=F)[,-3]
colnames(neuron)<-colnames(glia)<-c('gene','geneid')
neuron<-merge(FPKM,neuron,by.x='geneid',by.y='geneid')
glia<-merge(FPKM,glia,by.x='geneid',by.y='geneid')

write.table(file="neuron_TF_gene_expression.txt",neuron,sep="\t",quote=F,row.names=F)
write.table(file="glia_TF_gene_expression.txt",glia,sep="\t",quote=F,row.names=F)