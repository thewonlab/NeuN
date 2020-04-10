#"NeuN_Alzheimer_Module_Fishers"
options(stringsAsFactors=F)
library(dplyr)
library(ggplot2)

setwd("/proj/hyejunglab/NeuN/AD_H3K27ac/")
load("/proj/hyejunglab/chr/geneAnno_allgenes.rda")
geneAnno1<-geneAnno1[geneAnno1$hgnc_symbol!="" & geneAnno1$gene_biotype=="protein_coding" & geneAnno1$chromosome_name!="Y","hgnc_symbol"]
load("neg_hyper.genelist.rda")
load("neg_hypo.genelist.rda")
load("pos_hyper.genelist.rda")
load("pos_hypo.genelist.rda")
Module =read.csv("cels_208_mmc9.txt",header=T,sep=",")[,3:4]
module_list<-paste0("T-M", 1:20)
Module_gene=vector(length=20,mode="list")
for(i in c(1:20)){
    Module_gene[[i]]<-Module[Module[,2]==module_list[i],1]   #pay attention to this code!!!!
}
inputs<-list(neg_hyper.genelist, pos_hypo.genelist, neg_hypo.genelist, pos_hyper.genelist)#modification
id_sample<-c("neg","pos","neg","pos")#modification
ac_sample<-c("hyper","hypo","hypo","hyper")#modification
row_num<-length(geneAnno1)
metaMat<-matrix(0, nrow=row_num,ncol=3)
metaMat[,1]<-rownames(metaMat)<-geneAnno1  #assigned to gene symbol
metaMat<-data.frame(metaMat)

colnames(metaMat)<-c("reference","input_list","module_list")
fisher<-data.frame(matrix(ncol=5, nrow=length(Module_gene)*length(inputs)))#modification
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
        fisher[k,5]<-module_list[j]
        fisher[k,3]<- fisher.result$p.value
        fisher[k,4]<- fisher.result$estimate
        fisher[k,1]<- id_sample[i]#modification
        fisher[k,2]<- ac_sample[i]#modification
        metaMat[,2]<-metaMat[,3]<-c(0)
    }
}
#run fisher test for modules and hyper/hypo lists
colnames(fisher)<-c("pos_neg","hypo_hyper","p", "or",'Module')#modification
neg_hyper<-fisher[fisher[,1]==id_sample[1] & fisher[,2]==ac_sample[1],]
pos_hypo<-fisher[fisher[,1]==id_sample[2] & fisher[,2]==ac_sample[2],]
neg_hypo<-fisher[fisher[,1]==id_sample[3] & fisher[,2]==ac_sample[3],]
pos_hyper<-fisher[fisher[,1]==id_sample[4] & fisher[,2]==ac_sample[4],]

neg_hyper$FDR <- p.adjust(neg_hyper[,3],"fdr")
pos_hypo$FDR <- p.adjust(pos_hypo[,3],"fdr")
neg_hypo$FDR <- p.adjust(neg_hypo[,3],"fdr")
pos_hyper$FDR <- p.adjust(pos_hyper[,3],"fdr")

hyper<-rbind(neg_hyper,pos_hyper)
hypo<-rbind(neg_hypo,pos_hypo)

#r plot
# COLORS USED FOR EACH PLOT:
# neg_hyper - darkseagreen4
# neg_hypo - navy
# pos_hyper - darkseagreen2
# pos_hypo - lightblue

##########plot for all Modules
gplot<-function(input,sampleid,colorid){
    pdf(paste0("Modules_posneg_",sampleid,"_padj.pdf"), height=5, width=8)
    p<-ggplot(input, aes(x=Module,y=-log10(FDR), fill=pos_neg)) + theme_minimal() +coord_flip() + 
    geom_bar(stat="identity",position ="dodge") + scale_fill_manual(values=colorid)+ geom_hline(yintercept =-log10(0.01),color = "red")+
    labs(title=paste0(sampleid,"acetylated p-value (fdr corrected)"), y="-log10(FDR)", x="") + 
    theme_bw() + theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank())
    print(p)
    dev.off()
}
gplot(hypo,"hypo",c("navy", "lightblue"))
gplot(hyper,"hyper",c("darkseagreen4", "darkseagreen2"))

# Hypo for NeuNpos/neg
target_module<-c("T-M1", "T-M5", "T-M16", "T-M9", "T-M4","T-M7", "T-M8", "T-M18", "T-M10","T-M3", "T-M14")
plot2<-function(input,sampleid,colorid){
    output<-input %>% dplyr::filter(.,Module %in% target_module)
    gplot(output,sampleid,colorid)
}
plot2(hypo,"hypo_main",c("navy", "lightblue"))
plot2(hyper,"hyper_main",c("darkseagreen4", "darkseagreen2"))

##########################
#r plots for main figure
# Only NeuNpos-hypo, NeuNneg-hyper
hyper1 <- dplyr::filter(hyper, pos_neg=="neg") %>% 
        dplyr::mutate(., pos_neg.AD="NeuN- Hyperacetylated")

hypo1 <- dplyr::filter(hypo, pos_neg=="pos") %>%
        dplyr::mutate(., pos_neg.AD="NeuN+ Hypoacetylated")

joined <- dplyr::bind_rows(hyper1, hypo1)

joined$Module <-factor(joined$Module,levels =target_module)#Telling ggplot what order to put x-axis
gplot(joined,"hypo_hyper_all",c("navy", "lightblue"))

p2<-joined %>% dplyr::filter(.,Module %in% target_module)
gplot(p2,"hypo_hyper_filter",c("navy", "lightblue"))


