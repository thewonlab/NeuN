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

setwd("/nas/longleaf/home/hubenxia/project/Loops/neuron_vs_glia/")
load("geneAnno_allgenes.rda")
promoter <- read.table("gencode19_promoter.bed",sep="\t")
promoteranges <- GRanges(seqnames=promoter[,1], IRanges(promoter[,2],promoter[,3]), gene=promoter[,6])
NeuNneg_TSS <- read.table("NeuNneg_TSS_interacting_region_fixed.bed", header=T,sep="\t")
NeuNpos_TSS <- read.table("NeuNpos_TSS_interacting_region_fixed.bed", header=T,sep="\t")

tss <- function(input, promoteranges){
    result <- GRanges(seqnames=input[,1], IRanges(input[,2], input[,3]))
    mcols(result) <- input[ ,c(4:5)]
    olap <- findOverlaps(result, promoteranges)
    TSS_pro <- result[queryHits(olap)]
    mcols(TSS_pro) <- cbind(mcols(TSS_pro), mcols(promoteranges[subjectHits(olap)]))
    return(TSS_pro)
}

TSS_pro.NeuNpos <- tss(NeuNpos_TSS, promoteranges)
TSS_pro.NeuNneg <- tss(NeuNneg_TSS, promoteranges)

###############for overlapping interactions with differential H3K27ac peaks
neuron <- read.table("neuron_H3K27ac_differential_peak.txt",sep="\t",header=T,quote="")
glia <- read.table("non-neuron_H3K27ac_differential_peak.txt",sep="\t",header=T,quote="")

enhancer <- function(input, peak){
    result <- GRanges(seqnames=seqnames(input), IRanges(input$start_interaction, input$end_interaction))
    mcols(result) <- input$gene
    peak <- GRanges(seqnames=peak$Chrom, IRanges(peak$Start, peak$End))
    olap <- findOverlaps(peak,result)
    output <- peak[queryHits(olap)]
    mcols(output) <- mcols(result[subjectHits(olap)])
    return(unique(output))
}

pos_peak <- enhancer(TSS_pro.NeuNpos,neuron)
neg_peak <- enhancer(TSS_pro.NeuNneg,glia)

setwd("/proj/hyejunglab/epigenetics/NeuNsort/Panos/Loops_ATAC")
neuron_peak<-read.table("/proj/hyejunglab/epigenetics/NeuNsort/Panos/ATAC_neuron.bed",sep="\t",header=F)
glia_peak<-read.table("/proj/hyejunglab/epigenetics/NeuNsort/Panos/ATAC_glia.bed",sep="\t",header=F)

granges<-function(input){
    result<-GRanges(seqnames=input[,1], IRanges(input[,2],input[,3]))
    return(result)
}

neuronranges_peak<-granges(neuron_peak)
gliaranges_peak<-granges(glia_peak)

overlapping<-function(input,target,filename){
    olap<-findOverlaps(input,target)
    result<-unique(input[queryHits(olap)])
    print(head(result))
    write.table(file=paste0(filename,'_loops_ATACpeak.txt'),result,sep="\t",quote=F,row.names=F) #loops/ATAC-seq overlapping
}
overlapping(neuronranges_peak,pos_peak,'neuron')
overlapping(gliaranges_peak,neg_peak,'glia')