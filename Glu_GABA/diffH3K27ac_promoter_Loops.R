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


setwd("/nas/longleaf/home/hubenxia/project/Loops/GLU_GABA")

### Processing the data
load("/nas/longleaf/home/hubenxia/project/Loops/neuron_vs_glia/geneAnno_allgenes.rda")
promoter<-read.table("/nas/longleaf/home/hubenxia/project/Loops/neuron_vs_glia/gencode19_promoter.bed",sep="\t")
loops<-read.table("/proj/hyejunglab/NeuN/Loop/NeuNpos_TSS_interacting_region_fixed.bed", header=T,sep="\t")
colnames(loops)[c(4,5)]<-c('interacting_start', 'interacting_end')
Nloops<-nrow(loops)
### make GRanges
promoteranges<-GRanges(seqnames=promoter[,1], IRanges(promoter[,2],promoter[,3]), gene=promoter[,6])
loopsGR<-GRanges(seqnames=loops[,1], IRanges(loops[,2], loops[,3]))
mcols(loopsGR)<-loops[ ,c(4,5)]

### intersection of loops with promoter

overlap<-function(loop,promoteranges){
    olap<-findOverlaps(loop, promoteranges)
    loop_pro<-loop[queryHits(olap)]
    mcols(loop_pro)<-cbind(mcols(loop_pro), mcols(promoteranges[subjectHits(olap)]))
    ### output gene symbols
    return(loop_pro)
}
loops_pro<-overlap(loopsGR,promoteranges)
Nloops_pro<-nrow(unique(data.frame(loops_pro)[,c(1:3)]))


### make GRanges for differential H3K27ac peaks
peak <- read.table("/proj/hyejunglab/NeuN/ChIPseq/diffbind/H3K27ac3_differential_peaks.txt",sep="\t",header=T)
GLU_peak <- peak[peak$Fold >= 1 & peak$FDR<=0.01,]
GABA_peak <- peak[peak$Fold <= -1 & peak$FDR<=0.01,]

GLU_peak<-GRanges(seqnames=GLU_peak[,1], IRanges(GLU_peak[,2],GLU_peak[,3]))
GABA_peak<-GRanges(seqnames=GABA_peak[,1], IRanges(GABA_peak[,2],GABA_peak[,3]))


Analyze<-function(peak,loops,geneAnno1,flag){
    fetal_enh<-GRanges(seqnames=seqnames(loops),IRanges(loops$interacting_start,loops$interacting_end))
    mcols(fetal_enh)<-data.frame(data.frame(ranges(loops))[,1:2],gene=loops$gene)
    olap<-findOverlaps(fetal_enh, peak)
    olap2<-findOverlaps(peak, fetal_enh)
    magma <- peak[queryHits(olap2)]
    mcols(magma)<-cbind(mcols(magma), mcols(fetal_enh[subjectHits(olap2)]))
    magma <- unique(as.data.frame(magma)[,c(1:3,8)])
    fetal_enh<-fetal_enh[queryHits(olap)]
    mcols(fetal_enh)<-cbind(mcols(fetal_enh), mcols(peak[subjectHits(olap)]))
    result<-unique(data.frame(fetal_enh)[,c(1,6,7,2,3,8)])
    result<-merge(result,geneAnno1,by.x='gene',by.y='ensembl_gene_id')[,c(2:7)]
    result<-result[result$hgnc_symbol!='',]
    #write.table(file=paste0(flag,"_promoter_enhancer_magma_peak.txt"),magma,quote=F,sep="\t",row.names=F,col.names=F)
    write.table(file=paste0(flag,"_promoter_enhancer_looped_peak.txt"),unique(peak[subjectHits(olap)]),quote=F,sep="\t",row.names=F,col.names=F)
    write.table(file=paste0(flag,"_promoter_enhancer_loops_peak.txt"),result,quote=F,sep="\t",row.names=F,col.names=F)
    Nfetal_enh<-nrow(unique(data.frame(fetal_enh)[,1:3]))
    ### output gene symbols
    enh_gene<-unique(fetal_enh$gene)
    enh_gene<-geneAnno1[geneAnno1$ensembl_gene_id %in% enh_gene, "hgnc_symbol"]
    enh_gene<-enh_gene[enh_gene!=""]
    write.table(file=paste0(flag,"_promoter_enhancer_loops_peak_gene.txt"),unique(enh_gene),quote=F,sep="\t",row.names=F,col.names=F)
    return(list(loops=Nfetal_enh,gene=length(enh_gene)))
}

GLU<-Analyze(GLU_peak,loops_pro,geneAnno1,'GLUdiff')
GABA<-Analyze(GABA_peak,loops_pro,geneAnno1,'GABAdiff')

### make plots
# Nloops      187674
# Nloops_pro      15995
# NGLU  GLU$loops       9302
# NGABA GABA$loops      5448
# NGLUgene     GLU$gene     6238
# NGABAgene  GABA$gene 4344
# NGLU_peak GLU_peak 45911
# NGABA_peak GABA_peak 32169


id<-rep(c('GLU','GABA'),each=4)
region<-rep(c('promoter','nopromoter','enhancer','noenhancer'),2)
GLU_value<-c(Nloops_pro,Nloops-Nloops_pro,GLU$loops,Nloops-GLU$loops)
GABA_value<-c(Nloops_pro,Nloops-Nloops_pro,GABA$loops,Nloops-GABA$loops)
value<-c(GLU_value,GABA_value)
result<-data.frame(id,region,value)


gplot<-function(result,flag){
    pdf(paste0("loops_",flag,"_GLUvsGABA_peak.pdf"), width=5, height=5)
    p<-ggplot(result,aes(x=id,y=value,fill=region))+
    geom_bar(position="fill", stat="identity")+
    labs(title=paste0(flag," proportion"),y="proportion", x = flag)+
    theme_bw() + theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank())
    print(p)
    dev.off()
}
gplot(result[c(1:2,5:6),],'promoter')
gplot(result[c(3:4,7:8),],'enhancer')



pdf("loops_pro_enh_GLUvsGABA_peak.pdf", width=5, height=5)
ggplot(result,aes(x=id,y=value,fill=region))+
geom_bar(position="fill", stat="identity")+
labs(title="fetal proportion",y="proportion", x = "Hi-C samples")+
theme_bw() + theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank())
dev.off()