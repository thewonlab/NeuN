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
########################all poromoters
fire<-read.table("/nas/longleaf/home/hubenxia/project/FIRE/neuronvsglia/all_superFIRE.txt",sep="\t",header=T)
firediff.NeuNpos = fire[fire$pos==1,]
firediff.NeuNpos_number<-nrow(firediff.NeuNpos)
firediff.NeuNneg = fire[fire$neg==1,]
firediff.NeuNneg_number<-nrow(firediff.NeuNneg)
firanges.NeuNpos = GRanges(seqnames=firediff.NeuNpos$chr, IRanges(firediff.NeuNpos$start, firediff.NeuNpos$end))
firanges.NeuNneg = GRanges(seqnames=firediff.NeuNneg$chr, IRanges(firediff.NeuNneg$start, firediff.NeuNneg$end))

neuron<-read.table("neuron_H3K27ac_differential_peak.txt",sep="\t",header=T,quote="")
glia<-read.table("non-neuron_H3K27ac_differential_peak.txt",sep="\t",header=T,quote="")

neuronranges = GRanges(seqnames=neuron$Chrom, IRanges(neuron$Start, neuron$End))
gliaranges = GRanges(seqnames=glia$Chrom, IRanges(glia$Start, glia$End))

### [2] Overlap for NeuNpos with differetial neuronal and glial H3K27ac peaks 
olap = findOverlaps(firanges.NeuNpos, neuronranges)
firepeak.NeuNpos_neuron = firanges.NeuNpos[queryHits(olap)]
mcols(firepeak.NeuNpos_neuron) = cbind(mcols(firepeak.NeuNpos_neuron), mcols(neuronranges[subjectHits(olap)]))
firepeak.NeuNpos_neuron_DF<-annoGR2DF(firepeak.NeuNpos_neuron)
firepeak.NeuNpos_neuron_DF_fire<-firepeak.NeuNpos_neuron_DF[,c(1:3)]
firepeak.NeuNpos_neuron_DF_fire_number<-nrow(unique(firepeak.NeuNpos_neuron_DF_fire))


olap = findOverlaps(firanges.NeuNpos, gliaranges)
firepeak.NeuNpos_glia = firanges.NeuNpos[queryHits(olap)]
mcols(firepeak.NeuNpos_glia) = cbind(mcols(firepeak.NeuNpos_glia), mcols(gliaranges[subjectHits(olap)]))
firepeak.NeuNpos_glia_DF<-annoGR2DF(firepeak.NeuNpos_glia)
firepeak.NeuNpos_glia_DF_fire<-firepeak.NeuNpos_glia_DF[,c(1:3)]
firepeak.NeuNpos_glia_DF_fire_number<-nrow(unique(firepeak.NeuNpos_glia_DF_fire))


NeuNpos_nonpeak_neuron_number<-firediff.NeuNpos_number-firepeak.NeuNpos_neuron_DF_fire_number
NeuNpos_nonpeak_glia_number<-firediff.NeuNpos_number-firepeak.NeuNpos_glia_DF_fire_number

regions<-rep(c("peak","nonpeak","peak","nonpeak"),c(firepeak.NeuNpos_neuron_DF_fire_number,NeuNpos_nonpeak_neuron_number,firepeak.NeuNpos_glia_DF_fire_number,NeuNpos_nonpeak_glia_number))
sample_names<-rep(c("neuron","neuron","glia","glia"),c(firepeak.NeuNpos_neuron_DF_fire_number,NeuNpos_nonpeak_neuron_number,firepeak.NeuNpos_glia_DF_fire_number,NeuNpos_nonpeak_glia_number))

fire_peak_pos<-data.frame(count=regions,type=sample_names)


#######################plot
pdf("geom_bar_NeuNpos_allsuperFIRE_differentialH3K27ac_peak.pdf", width=5, height=5)
ggplot(fire_peak_pos,aes(x=type,fill=count))+
geom_bar(position ="fill")+
labs(title='proportion of all NeuNpos"\n"superFIREs overlapping with peak',y="proportion", x = "all superFIRE in NeuNpos")+
theme_bw() + theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank())
dev.off()




### [2] Overlap for NeuNneg with differetial neuronal and glial H3K27ac peaks 
olap = findOverlaps(firanges.NeuNneg, neuronranges)
firepeak.NeuNneg_neuron = firanges.NeuNneg[queryHits(olap)]
mcols(firepeak.NeuNneg_neuron) = cbind(mcols(firepeak.NeuNneg_neuron), mcols(neuronranges[subjectHits(olap)]))
firepeak.NeuNneg_neuron_DF<-annoGR2DF(firepeak.NeuNneg_neuron)
firepeak.NeuNneg_neuron_DF_fire<-firepeak.NeuNneg_neuron_DF[,c(1:3)]
firepeak.NeuNneg_neuron_DF_fire_number<-nrow(unique(firepeak.NeuNneg_neuron_DF_fire))


olap = findOverlaps(firanges.NeuNneg, gliaranges)
firepeak.NeuNneg_glia = firanges.NeuNneg[queryHits(olap)]
mcols(firepeak.NeuNneg_glia) = cbind(mcols(firepeak.NeuNneg_glia), mcols(gliaranges[subjectHits(olap)]))
firepeak.NeuNneg_glia_DF<-annoGR2DF(firepeak.NeuNneg_glia)
firepeak.NeuNneg_glia_DF_fire<-firepeak.NeuNneg_glia_DF[,c(1:3)]
firepeak.NeuNneg_glia_DF_fire_number<-nrow(unique(firepeak.NeuNneg_glia_DF_fire))


NeuNneg_nonpeak_neuron_number<-firediff.NeuNneg_number-firepeak.NeuNneg_neuron_DF_fire_number
NeuNneg_nonpeak_glia_number<-firediff.NeuNneg_number-firepeak.NeuNneg_glia_DF_fire_number

regions<-rep(c("peak","nonpeak","peak","nonpeak"),c(firepeak.NeuNneg_neuron_DF_fire_number,NeuNneg_nonpeak_neuron_number,firepeak.NeuNneg_glia_DF_fire_number,NeuNneg_nonpeak_glia_number))
sample_names<-rep(c("neuron","neuron","glia","glia"),c(firepeak.NeuNneg_neuron_DF_fire_number,NeuNneg_nonpeak_neuron_number,firepeak.NeuNneg_glia_DF_fire_number,NeuNneg_nonpeak_glia_number))

fire_peak_neg<-data.frame(count=regions,type=sample_names)


#######################plot
pdf("geom_bar_NeuNneg_all_superFIRE_differentialH3K27ac_peak.pdf", width=5, height=5)
ggplot(fire_peak_neg,aes(x=type,fill=count))+
geom_bar(position ="fill")+
labs(title='proportion of all NeuNneg"\n"superFIREs overlapping with peak',y="proportion", x = "all superFIRE in NeuNneg")+
theme_bw() + theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank())
dev.off()