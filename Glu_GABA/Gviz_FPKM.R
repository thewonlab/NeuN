options(stringsAsFactors=F)
library(ggplot2)
library(GenomicRanges)
library(plyr)
setwd("/nas/longleaf/home/hubenxia/project/Loops/GLU_GABA/")

FPKM<-unique(read.table("Glu_GABA_TPM.txt",header=T,sep="\t"))
target<-c("DCLK3","KCTD13","STAB1","GRIN2A","PTPRU","ATP2B2","SYNE1",'GRIK4','GAD1','ZNF804A','ITSN1','TTBK2','CNTN2','ADCY2')
FPKM<-FPKM[FPKM[,1] %in% target,]
rownames(FPKM)<-FPKM[,1]
FPKM<-FPKM[,-1]

FPKM_GABA<-FPKM[,c(1,3,5,7,9,11,13,15,17)]
FPKM_GLU<-FPKM[,c(2,4,6,8,10,12,14,16,18)]
FPKM<-cbind(FPKM_GABA,FPKM_GLU)
gplot<-function(input,target,filename){
    expr<-as.vector(t(input[rownames(input) %in% target,]))
    print(wilcox.test(expr[1:9],expr[10:18]))
    cell<-rep(c('GABA','Glu'),each=9)
    print(expr)
    expr<-data.frame(value=expr,celltype=cell)
    pdf(paste0(filename,"_Expression.pdf"))
    p<-ggplot(expr,aes(x=celltype, y=log(value), fill=celltype)) + 
        ylab(paste0('log(',target," expression)")) +  geom_boxplot(outlier.size = NA)  +
        ggtitle(paste0(filename," Expression"))+ theme_classic()
    print(p)
    dev.off()
}
gplot(FPKM,'GRIK4','GRIK4_TPM') 
gplot(FPKM,'GAD1','GAD1_TPM')
