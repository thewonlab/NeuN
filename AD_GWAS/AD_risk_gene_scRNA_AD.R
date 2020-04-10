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

option_list<-list(
    make_option(c("-a", "--anno"), type="character", default=NULL, 
            help="anno file name", metavar="character"),
    make_option(c("-n", "--neuron"), type="character", default=NULL, 
            help="neuron file name", metavar="character"),
    make_option(c("-g", "--glia"), type="character", default=NULL, 
            help="glia file name", metavar="character"),
    make_option(c("-o", "--output"), type="character", default=NULL, 
            help="output file name", metavar="character"),
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
if (is.null(opt$output)){
    print_help(opt_parser)
    stop("At least one argument must be supplied (output file).n", call.=FALSE)
}
if (is.null(opt$module)){
    print_help(opt_parser)
    stop("At least one argument must be supplied (module name).n", call.=FALSE)
}

load(opt$anno)
geneAnno1<-geneAnno1[geneAnno1$hgnc_symbol!="" & geneAnno1$gene_biotype=="protein_coding" & geneAnno1$chromosome_name!="Y","hgnc_symbol"]
#geneAnno1<-geneAnno1[geneAnno1$hgnc_symbol!="" & geneAnno1$chromosome_name!="Y","hgnc_symbol"]
neuron<-read.table(opt$neuron,sep="\t",header=F)[,1]
neuron<-neuron[neuron %in% geneAnno1]
glia<-read.table(opt$glia,sep="\t",header=F)[,1]
glia<-glia[glia %in% geneAnno1]

celldir = opt$module
filename = dir(celldir)  #list all files in this directory
filename = setdiff(filename,"genemodules.csv")  #remove
cellDEG = vector(length=length(filename),mode="list")

#################choose the DEGs.Ind.Model
for(i in 1:length(cellDEG)){
    file_id<-read.csv(paste0(celldir, filename[i]),header=T,sep=",",na.strings =c("#NAME?","NA",""))
    #file_id<-read.csv(paste0(celldir, filename[i]),header=T,na.strings =c("#NAME?","NA"," "),sep="\t") #early_late
    file_id_true_up<-file_id[file_id[,8]=="TRUE" & as.numeric(file_id[,5])>0.25,1]  #modification
    file_id_true_down<-file_id[file_id[,8]=="TRUE" & as.numeric(file_id[,5])<(-0.25),1]#modification
    file_id_true_up<-file_id_true_up[!is.na(file_id_true_up)]
    file_id_true_down<-file_id_true_down[!is.na(file_id_true_down)]
    cellDEG[[i]] = list(file_id_true_up,file_id_true_down)#modification
}
names(cellDEG)<-filename
###################

DEG_up<-lapply(cellDEG, '[[', 1)#modification
DEG_down<-lapply(cellDEG, '[[', 2)#modification
DEG<-list(DEG_up, DEG_down)#modification
inputs<-list(neuron,glia)#modification
cell_sample<-unlist(lapply(strsplit(names(cellDEG), split="_"),'[[',1))  #'[['===extraction 
pathology<- unlist(lapply(strsplit(names(cellDEG), split="_"),'[[',3)) #'[['===extraction
id_sample<-c("pos","neg")#modification

row_num<-length(geneAnno1)
metaMat<-matrix(0, nrow=row_num,ncol=3)
metaMat[,1]<-rownames(metaMat)<-geneAnno1  #assigned to gene symbol

metaMat<-data.frame(metaMat)
colnames(metaMat)<-c("reference","input_list","DEG")
regulation<-c("up","down")#modification
fisher<-data.frame(matrix(ncol=6, nrow=length(DEG_up)*length(DEG)*length(inputs)))#modification
k=0

for(i in c(1:length(inputs))){
    for(j in c(1:length(DEG))){
        for(m in c(1:length(DEG_up))){#modification
        metaMat[,2]<-ifelse(is.na(match(metaMat$reference,inputs[[i]])),0,1)
        metaMat[,3]<-ifelse(is.na(match(metaMat$reference,DEG[[j]][[m]])),0,1)#modification
        a1<-sum(metaMat$input_list==1 & metaMat$DEG==1)
        a2<-sum(metaMat$input_list==1 & metaMat$DEG==0)
        a3<-sum(metaMat$input_list==0 & metaMat$DEG==1)
        a4<-sum(metaMat$input_list==0 & metaMat$DEG==0)
        fisher.result<-fisher.test(matrix(c(a1,a2,a3,a4),2,2))
        k=k+1
        fisher[k,1]<- cell_sample[m]#modification
        fisher[k,2]<- fisher.result$p.value
        fisher[k,3]<- fisher.result$estimate
        fisher[k,4]<- id_sample[i]#modification
        fisher[k,5]<- pathology[m]#modification
        fisher[k,6]<- regulation[j]#modification
        metaMat[,2]<-metaMat[,3]<-c(0)
        }
    }
}

#run fisher test for modules and hyper/hypo lists
colnames(fisher)<-c("celltype", "p", "or","pos_neg","pathology","change")#modification
fisher.early = fisher[fisher$pathology=="earlypathologyVSnonpathology.csv",]
fisher.late = fisher[fisher$pathology=="pathologyVSnonpathology.csv",]

fisher.early$FDR <- p.adjust(fisher.early[,2],"fdr")
fisher.late$FDR <- p.adjust(fisher.late[,2],"fdr")

#############choose the DEGs.Ind.Model
write.table(file=paste0("fisher_Ind_",opt$output,".txt"),fisher,quote=F,sep="\t",row.names=F)
write.table(file=paste0("fisher.early_Ind_",opt$output,".txt"),fisher.early,quote=F,sep="\t",row.names=F)
write.table(file=paste0("fisher.late_Ind_",opt$output,".txt"),fisher.late,quote=F,sep="\t",row.names=F)
total<-rbind(fisher.early,fisher.late)
##########plot for all Modules
gplot<-function(input,sampleid,regulation,filename){
    input<-input[input$pos_neg==sampleid & input$change==regulation,]
    pdf(paste0(filename,"ADenrichment_padj.pdf"), height=5, width=8)
    p<-ggplot(input, aes(x=celltype,y=-log10(FDR))) + theme_minimal() +coord_flip() + 
    geom_bar(stat="identity",position ="dodge") + geom_hline(yintercept =-log10(0.01),color = "red")+
    labs(title=paste0(filename," p-value (fdr corrected)"), y="-log10(FDR)", x="") + 
    theme_bw() + theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank())
    print(p)
    dev.off()
}
gplot(fisher.early,"pos",'up',"neuron_early_ind")
gplot(fisher.early,"pos",'down',"neuron_early_ind")
gplot(fisher.early,"neg",'up',"glia_early_ind")
gplot(fisher.early,"neg",'down',"glia_early_ind")

gplot(fisher.late,"pos",'up',"neuron_late_ind_up")
gplot(fisher.late,"pos",'down',"neuron_late_ind_down")
gplot(fisher.late,"neg",'up',"glia_late_ind_up")
gplot(fisher.late,"neg",'down',"glia_late_ind_down")

gplot2<-function(input,sampleid,filename){
    input<-input[input$pos_neg==sampleid,]
    pdf(paste0(filename,"ADenrichment_padj.pdf"), height=5, width=8)
    p<-ggplot(input, aes(x=celltype,y=-log10(FDR),fill=pathology)) + theme_minimal() +coord_flip() + 
    geom_bar(stat="identity",position ="dodge") + geom_hline(yintercept =-log10(0.01),color = "red")+
    labs(title=paste0(filename," p-value (fdr corrected)"), y="-log10(FDR)", x="") + 
    theme_bw() + theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank())
    print(p)
    dev.off()
}
gplot2(total,"neg","glia_eary_late_ind")
gplot2(total,"pos","neuron_eary_late_ind")