######## [[1]] Promoter
### [1] Processing the data
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("Repitools")
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
library(Repitools)
#setwd("/proj/hyejunglab/NeuN/FIRE/")
load("geneAnno_allgenes.rda") # the saved file is geneAnno1

fire<-read.table("/nas/longleaf/home/hubenxia/project/FIRE/neuronvsglia/all_superFIRE.txt",sep="\t",header=T)
firediff.NeuNpos = fire[fire$pos==1,]
firediff.NeuNpos_number<-nrow(firediff.NeuNpos)
firediff.NeuNneg = fire[fire$neg==1,]
firediff.NeuNneg_number<-nrow(firediff.NeuNneg)
firanges.NeuNpos = GRanges(seqnames=firediff.NeuNpos$chr, IRanges(firediff.NeuNpos$start, firediff.NeuNpos$end))
firanges.NeuNneg = GRanges(seqnames=firediff.NeuNneg$chr, IRanges(firediff.NeuNneg$start, firediff.NeuNneg$end))


### [2] NeuNpos and NeuNneg overlapping with all promoter regions
promoter = read.table("gencode19_promoter.bed")
promoteranges = GRanges(seqnames=promoter[,1], IRanges(promoter[,2],promoter[,3]), gene=promoter[,6])

olap = findOverlaps(firanges.NeuNpos, promoteranges)
firepro.NeuNpos = firanges.NeuNpos[queryHits(olap)]
mcols(firepro.NeuNpos) = cbind(mcols(firepro.NeuNpos), mcols(promoteranges[subjectHits(olap)]))
firepro.NeuNpos_DF<-annoGR2DF(firepro.NeuNpos)
firepro.NeuNpos_DF_fire<-firepro.NeuNpos_DF[,c(1:3)]
firepro.NeuNpos_DF_fire_number<-nrow(unique(firepro.NeuNpos_DF_fire))

olap = findOverlaps(firanges.NeuNneg, promoteranges)
firepro.NeuNneg = firanges.NeuNneg[queryHits(olap)]
mcols(firepro.NeuNneg) = cbind(mcols(firepro.NeuNneg), mcols(promoteranges[subjectHits(olap)]))
firepro.NeuNneg_DF<-annoGR2DF(firepro.NeuNneg)
firepro.NeuNneg_DF_fire<-firepro.NeuNneg_DF[,c(1:3)]
firepro.NeuNneg_DF_fire_number<-nrow(unique(firepro.NeuNneg_DF_fire))
NeuNpos_nonpromoter_number<-firediff.NeuNpos_number-firepro.NeuNpos_DF_fire_number
NeuNneg_nonpromoter_number<-firediff.NeuNneg_number-firepro.NeuNneg_DF_fire_number
regions<-rep(c("promoter","nonpromoter","promoter","nonpromoter"),c(firepro.NeuNpos_DF_fire_number,NeuNpos_nonpromoter_number,firepro.NeuNneg_DF_fire_number,NeuNneg_nonpromoter_number))
sample_names<-rep(c("NeuNpos","NeuNpos","NeuNneg","NeuNneg"),c(firepro.NeuNpos_DF_fire_number,NeuNpos_nonpromoter_number,firepro.NeuNneg_DF_fire_number,NeuNneg_nonpromoter_number))

fire_promoter<-data.frame(count=regions,type=sample_names)
#firediff.NeuNpos_number:3966
#firediff.NeuNneg_number:3967
#firepro.NeuNpos_DF_fire_number:2329
#firepro.NeuNneg_DF_fire_number:2692
#NeuNpos_nonpromoter_number:1637
#NeuNneg_nonpromoter_number:1275
#######################plot
pdf("geom_bar_NeuN_all_superFIRE_promoter.pdf", width=5, height=5)
ggplot(fire_promoter,aes(x=type,fill=count))+
geom_bar(position ="fill")+
labs(title="proportion of all superFIREs overlapping with all Promoters",y="proportion", x = "all superFIRE")+
theme_bw() + theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank())
dev.off()