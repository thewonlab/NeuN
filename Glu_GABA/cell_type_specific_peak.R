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

setwd("/proj/hyejunglab/NeuN/ChIPseq/macs2")
Glu_peak<-read.table("Glu_peak.bed",sep="\t")
GABA_peak<-read.table("GABA_peak.bed",sep="\t")

GABA_peak<-GRanges(seqnames=GABA_peak[,1], IRanges(GABA_peak[,2],GABA_peak[,3]))
Glu_peak<-GRanges(seqnames=Glu_peak[,1], IRanges(Glu_peak[,2],Glu_peak[,3]))

###generated cell-type specifi ATAC-seq peaks
olap<-findOverlaps(Glu_peak, GABA_peak)
Glu_peak<-unique(Glu_peak[-queryHits(olap)])
GABA_peak<-unique(GABA_peak[-subjectHits(olap)])
write.table(file="Glu_specific.bed",Glu_peak,sep="\t",quote=F,row.names=F,col.names=T)
write.table(file="GABA_specific.bed",GABA_peak,sep="\t",quote=F,row.names=F,col.names=T)

