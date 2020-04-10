rm(list=ls())  #remove all of objects
options(stringsAsFactors=F)
library(ggplot2)
library(WGCNA)
library(reshape2)
library(dynamicTreeCut)
library(fastcluster)

# https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/

setwd("/nas/longleaf/home/hubenxia/project/Loops/neuron_vs_glia")
load("/proj/hyejunglab/chr/geneAnno_ensg_entrez_hgnc.rda")
reference <- geneAnno1[geneAnno1$hgnc_symbol !="",3]
posgene<-unique(read.table("inter_peak_NeuNposgene.txt",header=T)[,1])
neggene<-unique(read.table("inter_peak_NeuNneggene.txt",header=T)[,1])
targetgene<-list(posgene,neggene)
targetname<-c("NeuN+", "NeuN-")
listMat<-read.table("/proj/hyejunglab/traitsenrichment/Networks/NeelModule.txt",header=FALSE)

overlap <- function(target,input,refer,targetname,modulename){
    row_num<-length(refer)
    metaMat<-matrix(0, nrow=row_num,ncol=3)
    metaMat[,1]<-rownames(metaMat)<-refer  #assigned to gene symbol
    metaMat<-data.frame(metaMat)
    colnames(metaMat)<-c("reference","target","module")
    metaMat[,2]<-ifelse(is.na(match(metaMat$reference,target)),0,1)
    metaMat[,3]<-ifelse(is.na(match(metaMat$reference,input)),0,1)
    a1<-sum(metaMat$target==1 & metaMat$module==1)
    a2<-sum(metaMat$target==1 & metaMat$module==0)
    a3<-sum(metaMat$target==0 & metaMat$module==1)
    a4<-sum(metaMat$target==0 & metaMat$module==0)
    result<-fisher.test(matrix(c(a1,a2,a3,a4),2,2))
    return(c(result$p.value,result$estimate,targetname,modulename))
}

fisher<-data.frame(matrix(ncol=4, nrow=36))
m<-0
for (i in 1:length(targetname)){
    for(j in 1:18){
        modulename<-paste0("M",j)
        modulemat<-listMat[listMat[,3]==modulename, 2] #assigned to gene symbol
        modulemat<-modulemat[!duplicated(modulemat) & !is.na(modulemat)] #remove duplications and NA
        output<-overlap(targetgene[[i]],modulemat,reference,targetname[i],modulename)
        m<-m+1
        fisher[m,1] <- output[1]
        fisher[m,2] <- output[2]
        fisher[m,3] <- output[3]
        fisher[m,4] <- output[4]
    }
}
colnames(fisher)<-c('p','OR','cell','module')
fisher <- fisher[fisher$module!='M7',]
fisher$FDR <- p.adjust(fisher[,1],"fdr")
fisher$Enrichment <- -log10(fisher$FDR)
fisher$OR <- signif(as.numeric(fisher$OR),2)
fisher$Enrichment <- ifelse(fisher$OR>0, fisher$Enrichment, fisher$Enrichment * (-1))
fisher$cell <- factor(fisher$cell, levels=rev(c("NeuN+", "NeuN-")))
fisher$module <- factor(fisher$module, levels=paste0("M",c(1:6,8:18)))
##########select targeted modules Geneset
target<-c('M2','M3','M8','M9','M11','M13','M14','M16','M17','M18')
idname<-c('transcription','transcription','cell projection and development','ubiquitin','cell cyle',
        'neuronal and synaptic process','immune response','neuronal and synaptic process',
        'neuronal and synaptic process','immune response')
fisher_sub <- fisher[fisher$module %in% target,]
pdf("Enrichments_NeuNposvsneg_all_Mmodule_fisher.pdf",height=4,width=3)
theme_set(theme_bw())
ggplot(fisher, aes(x=cell, y=module)) + 
  geom_point(stat='identity', aes(col=-log10(FDR), size=OR)) + xlab("") + ylab("") + 
  scale_color_gradient2(low="#053061", mid="white", high="darkred") + theme(axis.text.x=element_text(angle = 90, hjust=1)) 
dev.off()
pdf("Enrichments_NeuNposvsneg_sub_Mmodule_fisher.pdf",height=4,width=3)
theme_set(theme_bw())
ggplot(fisher_sub, aes(x=cell, y=module)) + 
  geom_point(stat='identity', aes(col=-log10(FDR), size=OR)) + xlab("") + ylab("") + 
  scale_color_gradient2(low="#053061", mid="white", high="darkred") + theme(axis.text.x=element_text(angle = 90, hjust=1)) 
dev.off()
write.table(file="development_module_fisher.txt",fisher_sub,sep="\t",quote=F,row.names=F)
#####################
syn<-read.table("/proj/hyejunglab/traitsenrichment/CellTypes/synapse.sets")
colnames(syn)<-c("entrez", "funcgroup")
syn$ensg<-geneAnno1[match(syn$entrez, geneAnno1$entrezgene), "hgnc_symbol"]
syn<-syn[!is.na(syn$ensg),]
syngroup<-read.table("/proj/hyejunglab/traitsenrichment/CellTypes/synapse_clusters.index.txt", sep="\t", header=T)

fisher<-data.frame(matrix(ncol=4, nrow=36))
m<-0
for (i in 1:length(targetname)){
    for(j in 1:18){
        modulename<-paste0("func_gs",j)
        modulemat<-syn[syn[,2]==modulename, 3] #assigned to gene symbol
        modulemat<-modulemat[!duplicated(modulemat) & !is.na(modulemat)] #remove duplications and NA
        output<-overlap(targetgene[[i]],modulemat,reference,targetname[i],modulename)
        m<-m+1
        fisher[m,1] <- output[1]
        fisher[m,2] <- output[2]
        fisher[m,3] <- output[3]
        fisher[m,4] <- output[4]
    }
}
colnames(fisher)<-c('p','OR','cell','syn')
fisher$FDR <- p.adjust(fisher[,1],"fdr")
fisher$Enrichment <- -log10(fisher$FDR)
fisher$OR <- signif(as.numeric(fisher$OR),2)
fisher$Enrichment <- ifelse(fisher$OR>0, fisher$Enrichment, fisher$Enrichment * (-1))
fisher$cell <- factor(fisher$cell, levels=rev(c("NeuN+", "NeuN-")))
fisher$syn <- factor(fisher$syn, levels=paste0("func_gs",c(1:18)))

pdf("Enrichments_NeuNposvsneg_all_Msyn_fisher.pdf",height=5,width=3)
theme_set(theme_bw())
ggplot(fisher, aes(x=cell, y=syn)) + 
  geom_point(stat='identity', aes(col=-log10(FDR), size=OR)) + xlab("") + ylab("") + 
  scale_color_gradient2(low="#053061", mid="white", high="darkred") + theme(axis.text.x=element_text(angle = 90, hjust=1)) 
dev.off()
write.table(file="synapse_module_fisher.txt",fisher,sep="\t",quote=F,row.names=F)
