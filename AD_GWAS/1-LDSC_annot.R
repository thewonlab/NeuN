##Script to make annotation file for LD score regression
library(GenomicRanges)
library(GenomicFeatures)
library(stats4)
library(BiocGenerics)
library(parallel)
library(S4Vectors)
library(IRanges)
library(GenomeInfoDb)
library(AnnotationDbi)
library(Biobase)
library(gwascat)
library(Homo.sapiens)
library(AnnotationDbi)
library(OrganismDbi)
library(GO.db)
library(org.Hs.eg.db)
library(rtracklayer)
library(liftOver)

options(stringsAsFactors=FALSE)

setwd("/proj/hyejunglab/crossdisorder/LDSC/partherit")
outputdir<-"/proj/hyejunglab/crossdisorder/LDSC/partherit/annot_NeuN/"
neuron<-read.table("/nas/longleaf/home/hubenxia/project/Alzheimer/neuron_vs_glia/neuron_H3K27ac_differential_peak.txt",sep="\t",header=T,quote="")
glia<-read.table("/nas/longleaf/home/hubenxia/project/Alzheimer/neuron_vs_glia/non-neuron_H3K27ac_differential_peak.txt",sep="\t",header=T,quote="")
neuronranges<-GRanges(seqnames=neuron$Chrom, IRanges(neuron$Start, neuron$End))
gliaranges<-GRanges(seqnames=glia$Chrom, IRanges(glia$Start, glia$End))

newseqlevelnames<-substring(seqlevels(neuronranges),4,nchar(seqlevels(neuronranges)));
neuronranges<-renameSeqlevels(neuronranges,newseqlevelnames);

newseqlevelnames<-substring(seqlevels(gliaranges),4,nchar(seqlevels(gliaranges)));
gliaranges<-renameSeqlevels(gliaranges,newseqlevelnames);

neunpeakranges<-list(neuronranges,gliaranges)
neunmeta<-c("NeuNpos","NeuNneg")
## Check MHC region 
for (i in 1:22) {
    ##Read in locations of SNPs from the annotation files
    snps.table<-read.table(gzfile(paste0("/proj/hyejunglab/program/ldsc/LDSC/baselineLD_v1.1/baselineLD.",i,".annot.gz")),header=T)
    snps<-GRanges(snps.table$CHR,IRanges(snps.table$BP,snps.table$BP),
                  SNP=snps.table$SNP,CM=snps.table$CM)
    
    indicator<-matrix(0,length(snps),2)
    for(n in 1:length(neunmeta)){
        celltype<-neunmeta[n]
        crange<-neunpeakranges[[n]]
        olap.cat<-findOverlaps(snps,crange) #overlap
        indicator[queryHits(olap.cat),n] = 1   #queryHits output the index
    }
    print(paste0("NeuN done: chr",i))
    snps.cat<-snps
    mcols(snps.cat)<-cbind(as.data.frame(mcols(snps)),indicator)
    snps.cat.dat<-data.frame(snps.cat) #convert to a data.frame ("seqnames","start","end","width","strand","SNP","CM","X1","X2")
    outframe <- snps.cat.dat[,c("seqnames","start","SNP","CM","X1","X2")]
    colnames(outframe)<-c("CHR","BP","SNP","CM","NeuNpos","NeuNneg")
    print("#########")
    print(head(outframe))
    gz1<-gzfile(paste0(outputdir,"/",i,".annot.gz"), "w")
    write.table(outframe,gz1,quote=FALSE,col.names=TRUE,row.names=FALSE,sep="\t")
    close(gz1)
    print(i)
}