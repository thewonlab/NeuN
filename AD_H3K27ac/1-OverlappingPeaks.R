options(stringsAsFactors = F)
library(dplyr)
library(GenomicRanges)
library(ggplot2)
library(Vennerable)
library(stats4)
library(BiocGenerics)
library(parallel)
library(IRanges)

## [1] Finding Hypo/Hyper in each celltype
# Load differential H3K27ac peak files

neuron<-read.table("neuron_H3K27ac_differential_peak.txt",sep="\t",header=T,quote="")
glia<-read.table("non-neuron_H3K27ac_differential_peak.txt",sep="\t",header=T,quote="")
neuronranges<-GRanges(seqnames=neuron$Chrom, IRanges(neuron$Start, neuron$End))
gliaranges<-GRanges(seqnames=glia$Chrom, IRanges(glia$Start, glia$End))

# Load hypo/hyper files
hypo1<-read.csv("AD_hypoacetylated.csv")
hyper1<-read.csv("AD_hyperacetylated.csv")
load("geneAnno_allgenes.rda")

##modify chr23 and make GRanges
Grange<-function(peak){
    H3K27ac<- peak %>%
        dplyr::filter(., CHR!="24") %>%
        #choose rows/cases where conditions are true
        dplyr::mutate(CHR=paste0("chr",replace(CHR, CHR=="23", "X")) )%>%
        #mutate() adds new variables and preserves existing ones
        makeGRangesFromDataFrame(., keep.extra.columns = T)  
        #takes a data-frame-like object as input and tries to 
        #automatically find the columns that describe genomic ranges.
    return(H3K27ac)
}

hypo<-Grange(hypo1)
hyper<-Grange(hyper1)

#select only certain ranges from a GRanges which overlap something else
pos_hypo <- subsetByOverlaps(hypo, neuronranges)
length(pos_hypo)
neg_hypo <- subsetByOverlaps(hypo, gliaranges)
length(neg_hypo)
pos_hyper <- subsetByOverlaps(hyper, neuronranges)
length(pos_hyper)
neg_hyper <- subsetByOverlaps(hyper, gliaranges)
length(neg_hyper)
I obtained:                         Alex obtained
> length(pos_hypo)                  > length(pos_hypo)
[1] 1840                            [1] 1878
> length(neg_hypo)                  > length(neg_hypo)
[1] 691                             [1] 390
> length(pos_hyper)                 > length(pos_hyper)
[1] 73                              [1] 79
> length(neg_hyper)                 > length(neg_hyper)
[1] 937                             [1] 910

# plot proportional overlap between hypo/hyper files and neunpos/neg H3K27ac peaks
result <- data.frame(acetylation=c("hypo", "hypo", "hyper", "hyper"), 
                     celltype=c("neg", "pos", "neg", "pos"), 
                     count=c(length(neg_hypo)/length(hypo), 
                             length(pos_hypo)/length(hypo),
                             length(neg_hyper)/length(hyper),
                             length(pos_hyper)/length(hyper)))

pdf("Proportional_overlap.pdf", height=5, width=4)
ggplot(result,aes(x=acetylation,y=count,fill=celltype))+
geom_bar(stat="identity",position ="dodge")+ylim(0,1) + scale_fill_manual(values=c("navy", "lightblue"))+
labs(title="Proportion Overlap Between Hyper/Hypo\n and NeuNpos/NeuNneg Peaks",y="Proportional Overlap", x = "Acetylation Status")+
theme_bw() + theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank())
dev.off()

pdf("Proportional_overlap.pdf", height=5, width=8)
ggplot(result, aes(x=acetylation,y=count, fill=celltype)) + 
  geom_bar(stat="identity", position="dodge")+
  theme_minimal() +
  labs(title="Proportion Overlap Between Hyper/Hypo and NeuNpos/NeuNneg Peaks",
       x="Acetylation Status", y="Proportional Overlap") +
  ylim(0,1) + 
  scale_fill_manual(values=c("navy", "lightblue"))
dev.off()

#overlap acetylation peaks to promoter
#posTSS and negTSS to enhancer-centric
#load("Link_ATACregion2gene_NeuNposvsneg.rda")

#enhancer<-function(TSS){
#    enh <-data.frame(TSS) %>% dplyr::mutate(seqnames=paste0("chr",seqnames)) %>%  
#     #convert GRanges to data.frame
#        dplyr::mutate(., enh_end= enhancer+10000, TSS=start) %>% 
#        dplyr::select(., seqnames, enhancer, enh_end, TSS, gene) %>%
#        dplyr::rename(.,start=enhancer, end=enh_end) %>%
#        makeGRangesFromDataFrame(., keep.extra.columns = T)
#    return(enh)
#}
#pos_enh<-enhancer(posTSSint)
#neg_enh<-enhancer(negTSSint)
#seqnames, enhancer, enh_end, TSS, gene(1,5,6,2,7)

TSS_pos<-unique(read.table("TSS_pro.NeuNpos_DF.txt",header=T,sep="\t"))
TSS_neg<-unique(read.table("TSS_pro.NeuNneg_DF.txt",header=T,sep="\t"))

pos_enh<-GRanges(seqnames=TSS_pos[,1], IRanges(start=TSS_pos[,5], end=TSS_pos[,6]),TSS=TSS_pos[,2],gene=TSS_pos[,7])
neg_enh<-GRanges(seqnames=TSS_neg[,1], IRanges(start=TSS_neg[,5], end=TSS_neg[,6]),TSS=TSS_neg[,2],gene=TSS_neg[,7])


# SubsetByOverlaps to find Hypo/Hyper in NeuNpos/neg
suboverlap<-function(enh,peak){
    enh_peak_overlap<-subsetByOverlaps(enh,peak) %>% data.frame(.)
    return(enh_peak_overlap)
}

neunpos_hypo<-suboverlap(pos_enh, pos_hypo)
neunpos_hyper<-suboverlap(pos_enh, pos_hyper)
neunneg_hypo<-suboverlap(neg_enh, neg_hypo)
neunneg_hyper<-suboverlap(neg_enh, neg_hyper)


save(neunpos_hypo, file="neunpos_hypo.rda")
save(neunpos_hyper, file="neunpos_hyper.rda")
save(neunneg_hypo, file="neunneg_hypo.rda")
save(neunneg_hyper, file="neunneg_hyper.rda")

#make genelists from enhancer-gene interactions
# Compile genelists
pos_hyper.genelist <- unique(neunpos_hyper$gene)
save(pos_hyper.genelist, file= "pos_hyper.genelist.rda")
neg_hyper.genelist <- unique(neunneg_hyper$gene)
save(neg_hyper.genelist, file="neg_hyper.genelist.rda")
pos_hypo.genelist <- unique(neunpos_hypo$gene)
save(pos_hypo.genelist, file="pos_hypo.genelist.rda")
neg_hypo.genelist <- unique(neunneg_hypo$gene)
save(neg_hypo.genelist, file="neg_hypo.genelist.rda")


#r Venn Diagrams

Vennplot<-function(pos_gene,neg_gene,sample_name){
    overlap_gene<-length(dplyr::intersect(pos_gene,neg_gene)) #find common entries
    print(overlap_gene)
    pos_spec<-nrow(pos_gene)-overlap_gene# NeuN+ specific
    neg_spec<-nrow(neg_gene)-overlap_gene# NeuN- specific
    venn<-list("NeuNpos"= pos_gene,"NeuNneg"=neg_gene)
    VD<-Venn(venn)
    pdf(paste0(sample_name,"_Venn_Diagram",".pdf"), height=5, width=5)
    plot(VD,doWeights=T)
    dev.off()
}

length(dplyr::intersect(pos_hypo.genelist, neg_hypo.genelist))

Vennplot(pos_hypo.genelist, neg_hypo.genelist,"Hypoacetylated")
Vennplot(pos_hyper.genelist, neg_hyper.genelist,"Hyperacetylated")