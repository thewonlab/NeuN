### [1] Processing the data
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

setwd("/nas/longleaf/home/hubenxia/project/Loops/neuron_vs_glia")
load("geneAnno_allgenes.rda") # the saved file is geneAnno1
promoter<-read.table("gencode19_promoter.bed",sep="\t")
glia<-read.table("NeuNneg_TSS_interacting_region_fixed.bed", header=T,sep="\t")
glia_number<-nrow(glia)
#glia_number:167551
neuron<-read.table("NeuNpos_TSS_interacting_region_fixed.bed", header=T,sep="\t")
neuron_number<-nrow(neuron)
#neuron_number:187674

promoteranges<-GRanges(seqnames=promoter[,1], IRanges(promoter[,2],promoter[,3]), gene=promoter[,6])

### [2] Overlap  for interactions with promoters

overlap<-function(cell,promoteranges,flag){
    TSS_loops<-GRanges(seqnames=cell[,1], IRanges(cell[,2], cell[,3]))
    mcols(TSS_loops)<-cell[ ,c(4:5)]
    olap<-findOverlaps(TSS_loops, promoteranges)
    TSS_pro<-TSS_loops[queryHits(olap)]
    TSS_pro_num<-nrow(unique(as.data.frame(TSS_pro)))
    TSS<-as.data.frame(mcols(TSS_pro))
    inter_pro<-GRanges(seqnames=seqnames(TSS_pro), IRanges(TSS_pro$start_interaction,TSS_pro$end_interaction))
    mcols(inter_pro)<-cbind(TSS, mcols(promoteranges[subjectHits(olap)]))

    olap<-findOverlaps(promoteranges,inter_pro)
    pro<-unique(as.data.frame(promoteranges[queryHits(olap)]))
    write.table(file=paste0(flag,"_promoter_promoter_seq.txt"),pro,sep="\t",quote=F,row.names=F,col.names=F)

    olap<-findOverlaps(inter_pro, promoteranges)
    inter_pro_pro<-inter_pro[queryHits(olap)]
    write.table(file=paste0(flag,"_promoter_promoter_interaction.txt"),as.data.frame(inter_pro_pro),sep="\t",quote=F,row.names=F)
    result<-unique(as.data.frame(inter_pro_pro)[,1:7])
    write.table(file=paste0(flag,"_promoter_promoter_interaction_unique.txt"),result,sep="\t",quote=F,row.names=F)
    return(list(TSS_pro_num-nrow(result),nrow(result)))
}
pro_pro_neuron<-overlap(neuron,promoteranges,'neruon')
pro_pro_glia<-overlap(glia,promoteranges,'glia')

result<-data.frame(c(pro_pro_glia[[1]],pro_pro_glia[[2]]),c(pro_pro_neuron[[1]],pro_pro_neuron[[2]]),c('other-p','p-p'))
colnames(result)<-c('glia','neuron','type')
result<-melt(result)
     type variable value
1 other-p     glia 60798
2     p-p     glia 19561
3 other-p   neuron 70460
4     p-p   neuron 21558
pdf("geom_bar_promoter_promoter.pdf", width=5, height=5)
ggplot(result,aes(x=factor(variable,levels=c('neuron','glia')),y=value,fill=type))+
    geom_bar(stat="identity",position = "fill")+
    labs(title="proportion of P-P \n other_p_pos:70460 \n opter_p_neg:60798 \n p_p_pos:21558 \n p_p_neg:19561",y="proportion")+
    theme_classic()
dev.off()

setwd("/nas/longleaf/home/hubenxia/project/Loops/neuron_vs_glia")
neruon<-read.table("neruon_promoter_promoter_seq.txt",sep="\t",header=T)
glia<-read.table("glia_promoter_promoter_seq.txt",sep="\t",header=T)
neuron_atac<-read.csv("/proj/hyejunglab/epigenetics/NeuNsort/Panos/NeuNpos_ATAC_Panos.csv")[,1:3]
glia_atac<-read.csv("/proj/hyejunglab/epigenetics/NeuNsort/Panos/NeuNneg_ATAC_Panos.csv")[,1:3]

setwd("/proj/hyejunglab/NeuN/Loop/p_p_loops")
overlapping<-function(input,target,filename){
    input<-GRanges(seqnames=input[,1], IRanges(input[,2],input[,3]))
    target<-GRanges(seqnames=target[,1], IRanges(target[,2],target[,3]))
    olap<-findOverlaps(input,target)
    result<-unique(input[queryHits(olap)])
    write.table(file=paste0(filename,'_ATAC_PP.txt'),result,sep="\t",quote=F,row.names=F)
}
overlapping(neuron_atac,neruon,'neuron')
overlapping(glia_atac,glia,'glia')
