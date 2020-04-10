#!/usr/bin/env Rscript
library(stats4)
library(BiocGenerics)
library(parallel)
library(S4Vectors)
library(reshape2)
library(IRanges)
library(GenomeInfoDb)
library(optparse)
options(stringsAsFactors=F)
library(dplyr)
library(ggplot2)
library(plyr)

#"NeuN_Alzheimer_Module_Fishers"
option_list<-list(
    make_option(c("-a", "--anno"), type="character", default=NULL, 
            help="anno file name", metavar="character"),
    make_option(c("-n", "--neuron"), type="character", default=NULL, 
            help="neuron file name", metavar="character"),
    make_option(c("-g", "--glia"), type="character", default=NULL, 
            help="glia file name", metavar="character"),
    make_option(c("-m", "--module"), type="character", default=NULL, 
            help="module file name", metavar="character")
)

opt_parser<-OptionParser(option_list=option_list)
opt<-parse_args(opt_parser)

if (is.null(opt$anno)){
    print_help(opt_parser)
    stop("At least one argument must be supplied (anno file).n", call.=FALSE)
}
if (is.null(opt$neuron)){
    print_help(opt_parser)
    stop("At least one argument must be supplied (neuron file).n", call.=FALSE)
}
if (is.null(opt$glia)){
    print_help(opt_parser)
    stop("At least one argument must be supplied (glia file).n", call.=FALSE)
}
if (is.null(opt$module)){
    print_help(opt_parser)
    stop("At least one argument must be supplied (module name).n", call.=FALSE)
}


load(opt$anno)
#geneAnno1<-geneAnno1[geneAnno1$hgnc_symbol!="" & geneAnno1$gene_biotype=="protein_coding" & geneAnno1$chromosome_name!="Y","hgnc_symbol"]
geneAnno1<-geneAnno1[geneAnno1$hgnc_symbol!="" & geneAnno1$chromosome_name!="Y","hgnc_symbol"]
neuron<-read.table(opt$neuron,sep="\t",header=F)[,1]
glia<-read.table(opt$glia,sep="\t",header=F)[,1]
Module =read.csv(opt$module,header=T,sep=",")[,3:4]
module_list<-paste0("T-M", 1:20)
Module_gene=vector(length=20,mode="list")
for(i in c(1:20)){
    Module_gene[[i]]<-Module[Module[,2]==module_list[i],1]   #pay attention to this code!!!!
}
inputs<-list(neuron,glia)#modification
id_sample<-c("pos","neg")#modification
row_num<-length(geneAnno1)
metaMat<-matrix(0, nrow=row_num,ncol=3)
metaMat[,1]<-rownames(metaMat)<-geneAnno1  #assigned to gene symbol
metaMat<-data.frame(metaMat)

colnames(metaMat)<-c("reference","input_list","module_list")
fisher<-data.frame(matrix(ncol=4, nrow=length(Module_gene)*length(inputs)))#modification
k=0
for(i in c(1:length(inputs))){
    for(j in c(1:length(Module_gene))){
        metaMat[,3]<-ifelse(is.na(match(metaMat$reference,inputs[[i]])),0,1)
        metaMat[,2]<-ifelse(is.na(match(metaMat$reference,Module_gene[[j]])),0,1)
        a1<-sum(metaMat$input_list==1 & metaMat$module_list==1)
        a2<-sum(metaMat$input_list==1 & metaMat$module_list==0)
        a3<-sum(metaMat$input_list==0 & metaMat$module_list==1)
        a4<-sum(metaMat$input_list==0 & metaMat$module_list==0)
        fisher.result<-fisher.test(matrix(c(a1,a2,a3,a4),2,2))
        k=k+1
        fisher[k,4]<-module_list[j]
        fisher[k,3]<- fisher.result$p.value
        fisher[k,2]<- fisher.result$estimate
        fisher[k,1]<- id_sample[i]#modification
        metaMat[,2]<-metaMat[,3]<-c(0)
    }
}
#run fisher test for modules
colnames(fisher)<-c("pos_neg","p", "or",'Module')#modification
pos<-fisher[fisher[,1]==id_sample[1],]
neg<-fisher[fisher[,1]==id_sample[2],]

pos$FDR <- p.adjust(pos[,3],"fdr")
neg$FDR <- p.adjust(neg[,3],"fdr")

pos_neg<-rbind(pos,neg)
write.table(file="neuron_AD_module.txt",pos,sep="\t",quote=F,row.names=F)
write.table(file="glia_AD_module.txt",neg,sep="\t",quote=F,row.names=F)
write.table(file="neuron_glia_AD_module.txt",pos_neg,sep="\t",quote=F,row.names=F)
##########plot for all Modules
gplot<-function(input,sampleid,colorid){
    pdf(paste0("Modules_",sampleid,"_padj.pdf"), height=5, width=8)
    p<-ggplot(input, aes(x=Module,y=-log10(FDR), fill=pos_neg)) + theme_minimal() +coord_flip() + 
    geom_bar(stat="identity",position ="dodge") + scale_fill_manual(values=colorid)+ geom_hline(yintercept =-log10(0.01),color = "red")+
    labs(title=paste0(sampleid," p-value (fdr corrected)"), y="-log10(FDR)", x="") + 
    theme_bw() + theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank())
    print(p)
    dev.off()
}
gplot(pos,"posAD",c("navy", "lightblue"))
gplot(neg,"negAD",c("darkseagreen4", "darkseagreen2"))

# Hypo for NeuNpos/neg
target_module<-c("T-M1", "T-M5", "T-M16", "T-M9", "T-M4","T-M7", "T-M8", "T-M18", "T-M10","T-M3", "T-M14")
plot2<-function(input,sampleid,target_module,colorid){
    output<-input %>% dplyr::filter(.,Module %in% target_module)
    output$Module <-factor(output$Module,levels =target_module)
    gplot(output,sampleid,colorid)
}

plot2(pos,"posAD_main",target_module,c("navy", "lightblue"))
plot2(neg,"negAD_main",target_module,c("darkseagreen4", "darkseagreen2"))

##########################
#r plots for main figure

hyper1 <- dplyr::filter(pos, pos_neg=="pos") %>% 
        dplyr::mutate(., pos_neg.AD="posAD")

hypo1 <- dplyr::filter(neg, pos_neg=="neg") %>%
        dplyr::mutate(., pos_neg.AD="negAD")

joined <- dplyr::bind_rows(hyper1, hypo1)

joined$Module <-factor(joined$Module,levels =target_module) #Telling ggplot what order to put x-axis
gplot(joined,"pos_neg_AD_all",c("navy", "lightblue"))

p2<-joined %>% dplyr::filter(.,Module %in% target_module)
gplot(p2,"pos_neg_AD_filter",c("navy", "lightblue"))