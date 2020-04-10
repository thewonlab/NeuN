#module add r/3.6.0
#R
library(DiffBind)
library(GenomicRanges)
library(stats4)
library(BiocGenerics)
library(parallel)
library(S4Vectors)
library(IRanges)
library(GenomeInfoDb)
library(SummarizedExperiment)
library(Biobase)
library(pheatmap)
library(DelayedArray)
library(BiocParallel)
library(matrixStats)
options(stringsAsFactors=FALSE)

fnamesamplesheet <- "/nas/longleaf/home/hubenxia/project/ChIP-seq/GLU_GABA/GLU_GABA.csv"
blacklist <- read.table("/proj/steinlab/projects/atacseq/MinimalNecessary_Dan/BlackListhg19/wgEncodeHg19ConsensusSignalArtifactRegions.bed",header=FALSE,sep="\t")
colnames(blacklist)[1:3] <- c("chrom","chromStart","chromEnd")
blacklist <- GRanges(blacklist$chrom, IRanges(blacklist$chromStart,blacklist$chromEnd))

##Remove all peaks which map to ENCODE blacklisted regions
samplesheet <- read.csv(fnamesamplesheet,header=TRUE)
for (i in 1:length(samplesheet$Peaks)) {
  broadPeak <- read.table(samplesheet$Peaks[i],header=FALSE,sep="\t")
  colnames(broadPeak) <- c("chrom","chromStart","chromEnd","name","score","strand +/-","signalValue","pValue","qValue")
  #broadPeak$chrom <- unlist(strsplit(broadPeak$chrom, split="chr"))[seq(2,2*length(broadPeak$chrom),2)]
  #broadPeak <- broadPeak[broadPeak$chrom %in% c(paste0('chr',1:22),"chrX"),]
  broadPeakbed <- GRanges(broadPeak$chrom,IRanges(broadPeak$chromStart,broadPeak$chromEnd))
  mcols(broadPeakbed) <- broadPeak[,4:9]
  olap <- findOverlaps(broadPeakbed,blacklist)
  broadPeakbedfilt <- broadPeakbed[-queryHits(olap)]
  ##Write out the filtered peaks
  broadPeakbedfiltout <- data.frame(chrom=seqnames(broadPeakbedfilt),chromStart=start(broadPeakbedfilt),chromEnd=end(broadPeakbedfilt))
  broadPeakbedfiltout <- cbind(broadPeakbedfiltout,as.data.frame(mcols(broadPeakbedfilt)));
  ##Make a new output file name
  pos <- unlist(strsplit(basename(samplesheet$Peaks[i]),split='.',fixed=T))[1:3]
  filename=paste(pos[1],pos[2],pos[3],"noblacklist.bed",sep='.')
  write.table(broadPeakbedfiltout,file=paste0("/proj/hyejunglab/NeuN/ChIPseq/cleanpeaks/",filename),sep="\t",col.names=FALSE,row.names=FALSE,quote=FALSE)
}


setwd("/proj/hyejunglab/NeuN/ChIPseq/diffbind")
input<- dba(sampleSheet="/nas/longleaf/home/hubenxia/project/ChIP-seq/GLU_GABA/GLU_GABA_clean.csv",bCorPlot=FALSE,minOverlap=0.4)
pdf("Correlation_heatmap_occupancy_peak_caller_score.pdf")
plot(input)
dev.off()

input<-dba.count(input, summits=0)
pdf("Correlation_heatmap_affinity_read_count.pdf")
plot(input)
dev.off()
pdf("PCA_plot_affinity_for_all_sites.pdf")
dba.plotPCA(input,DBA_TISSUE,label=DBA_CONDITION)
dev.off()

#input<- dba.contrast(input, categories=DBA_CONDITION)
input <- dba.contrast(input, group1=input$masks$Glu, group2=input$masks$GABA, name1="Glu",name2="GABA", block=DBA_REPLICATE)

#snames <- c(colnames(input$class[,input$masks$Glu]),colnames(input$class[,input$masks$GABA]))
#(length(unique(snames))!=sum(input$masks$Glu|input$masks$GABA))

input <- dba.analyze(input)
pdf("Correlation_heatmap_significantly_differentially_bound_sites.pdf")
plot(input, contrast=1)
dev.off()
pdf("PCA_plot_affinity_differentially_bound_sites.pdf")
dba.plotPCA(input,contrast=1,label=DBA_TISSUE)
dev.off()
pdf("MA_plot_GLU_GABA_contrast.pdf")
dba.plotMA(input)
dev.off()
pdf("Volcano_plot_GLU_GABA_contrast.pdf")
dba.plotVolcano(input)
dev.off()

input_differential_bind<- dba.report(input)
print(head(input_differential_bind))
write.table(input_differential_bind,file="H3K27ac3_differential_peaks.txt",row.names=F,quote = F,sep="\t")
