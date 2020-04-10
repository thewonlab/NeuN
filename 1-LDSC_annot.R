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
options(stringsAsFactors=FALSE)

setwd("/proj/hyejunglab/crossdisorder/LDSC/partherit")
outputdir<-"/proj/hyejunglab/crossdisorder/LDSC/partherit/annot_Glu_GABA/"
peak <- read.table("/proj/hyejunglab/NeuN/ChIPseq/diffbind/H3K27ac3_differential_peaks.txt",sep="\t",header=T)
GLU_peak <- peak[peak$Fold >= 1 & peak$FDR<=0.01,][,c(1:3)]
GABA_peak <- peak[peak$Fold <= -1 & peak$FDR<=0.01,][,c(1:3)]

granges<-function(input_file){
    input_ranges<-GRanges(seqnames=input_file[,1], IRanges(input_file[,2],input_file[,3]))
    return(input_ranges)
}
GLUranges<-granges(GLU_peak)
GABAranges<-granges(GABA_peak)

rm_chr<-function(input){
    newseqlevelnames<-substring(seqlevels(input),4,nchar(seqlevels(input)));
    input<-renameSeqlevels(input,newseqlevelnames)
    return(input)
}

GLUranges<-rm_chr(GLUranges)
GABAranges<-rm_chr(GABAranges)

ATAC<-list(GLUranges,GABAranges)
neunmeta<-c("GLU","GABA")
## Check MHC region 
for (i in 1:22) {
    ##Read in locations of SNPs from the annotation files
    snps.table<-read.table(gzfile(paste0("/proj/hyejunglab/program/ldsc/LDSC/baselineLD_v1.1/baselineLD.",i,".annot.gz")),header=T)
    snps<-GRanges(snps.table$CHR,IRanges(snps.table$BP,snps.table$BP),
                  SNP=snps.table$SNP,CM=snps.table$CM)
    indicator<-matrix(0,length(snps),2)
    for(n in 1:length(neunmeta)){
        celltype<-neunmeta[n]
        crange<-ATAC[[n]]
        olap.cat<-findOverlaps(snps,crange) #overlap
        indicator[queryHits(olap.cat),n] = 1   #queryHits output the index
    }
    print(paste0("NeuN done: chr",i))
    snps.cat<-snps
    mcols(snps.cat)<-cbind(as.data.frame(mcols(snps)),indicator)
    snps.cat.dat<-data.frame(snps.cat) #convert to a data.frame ("seqnames","start","end","width","strand","SNP","CM","X1","X2")
    outframe <- snps.cat.dat[,c("seqnames","start","SNP","CM","X1","X2")]
    colnames(outframe)<-c("CHR","BP","SNP","CM","GLU","GABA")
    gz1<-gzfile(paste0(outputdir,"/",i,".annot.gz"), "w")
    write.table(outframe,gz1,quote=FALSE,col.names=TRUE,row.names=FALSE,sep="\t")
    close(gz1)
    print(i)
}