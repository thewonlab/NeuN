options(stringsAsFactors=F)
library(ggplot2)
library(GenomicRanges)
library(plyr)
setwd("/nas/longleaf/home/hubenxia/project/RNA-seq/neuron_vs_glia/known_gene/")
load("/proj/hyejunglab/NeuN/RNAseq/analysis/datExprMeta.rda")
FPKM<-read.table("neuron_glia_FPKM_dot_remove.txt",header=T,sep="\t",row.names=1)
colnames(FPKM)<-c("CMC_HBCC_028_neg","CMC_HBCC_025_pos","CMC_HBCC_028_pos","CMC_HBCC_027_pos",
                        "CMC_HBCC_025_neg","CMC_HBCC_027_neg","CMC_HBCC_016_pos","CMC_HBCC_016_neg")

FPKM_neuron<-FPKM[,c(2:4,7)]
FPKM_glia<-FPKM[,c(1,5:6,8)]
FPKM<-cbind(FPKM_neuron,FPKM_glia)
EHD1<-'ENSG00000110047'
CACNG3<-'ENSG00000006116'

gplot<-function(input,target,flag,filename){
    expr<-as.vector(t(input[rownames(input) %in% target,]))
    print(expr)
    print(wilcox.test(expr[1:4],expr[5:8]))
    cell<-rep(c('NeuN+','NeuN-'),each=4)
    expr<-data.frame(value=expr,celltype=cell)
    pdf(paste0(filename,"_Expression.pdf"))
    p<-ggplot(expr,aes(x=celltype, y=log(value), fill=celltype)) + 
        ylab(paste0('log(',flag," expression)")) +  geom_boxplot(outlier.size = NA)  +
        ggtitle(paste0(filename," Expression"))+ theme_classic()
    print(p)
    dev.off()
}
gplot(FPKM,EHD1,'EHD1','EHD1_FPKM') #
gplot(FPKM,CACNG3,'CACNG3','CACNG3_FPKM')