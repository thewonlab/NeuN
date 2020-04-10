rm(list=ls())
options(stringsAsFactors=F)
library(ggplot2)
library(reshape2)
library(scales)
setwd("/nas/longleaf/home/hubenxia/project/ATAC-seq")
neuron<-read.table("neuron_FPKM_TF.txt",sep="\t",header=T)
glia<-read.table("glia_FPKM_TF.txt",sep="\t",header=T)
Glia_neuron<-glia[,c(2,3,5)]
colnames(Glia_neuron)<-c('gene','FPKM','score')
Glia_neuron$cell<-'neuron'

Glia_glia<-glia[,c(2,4,6)]
colnames(Glia_glia)<-c('gene','FPKM','score')
Glia_glia$cell<-'glia'

Neuron_neuron<-neuron[,c(2,3,5)]
colnames(Neuron_neuron)<-c('gene','FPKM','score')
Neuron_neuron$cell<-'neuron'
Neuron_glia<-neuron[,c(2,4,6)]
colnames(Neuron_glia)<-c('gene','FPKM','score')
Neuron_glia$cell<-'glia'

motif<-Reduce(function(x, y) rbind(x, y), list(Glia_neuron,Glia_glia,Neuron_neuron,Neuron_glia))
#motif<-read.table("neuron_glia_motif.txt",sep="\t",header=T)
motif$cell <- factor(motif$cell, levels=c("neuron", "glia"))
mid<- 0

id<-c(neuron[,2],glia[,2])
pdf('motif_expression.pdf',width=3,height=6)
theme_set(theme_bw())
p1<-ggplot(motif, aes(x=factor(gene,levels=id), y=cell))+ coord_flip() + 
    geom_point(stat='identity', aes(size=score, col=log(FPKM))) + xlab("") + ylab("") + 
    theme(axis.text.x=element_text(angle = 90, hjust=1))+
    scale_color_gradient2(low="lightblue", mid="white", high="darkred")
print(p1)
dev.off()


pdf('motif_expression2.pdf',width=2.5,height=4)
theme_set(theme_bw())
ggplot(motif, aes(x=factor(gene,levels=id), y=cell))+ coord_flip() + 
    geom_point(stat='identity', aes(col=score, size=log(FPKM))) + xlab("") + ylab("") + 
    theme(axis.text.x=element_text(angle = 90, hjust=1))+
    scale_color_gradient2(low="lightblue", mid="white", high="darkred")
dev.off()