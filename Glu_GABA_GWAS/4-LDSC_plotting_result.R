options(stringsAsFactors=F)
library(WGCNA)
library(reshape)
library(ggplot2)
setwd("/proj/hyejunglab/crossdisorder/LDSC/partherit/results/phase3/Glu_GABA/")
directory<-list.files(pattern = ".*\\.results")
conditions<-unlist(strsplit(directory, split="_results[.]results"))

condition<-conditions[1]
results<-read.table(paste0(condition,"_results.results"), header=T)
grepneun<-grep("^GLU|^GABA", results$Category)

neun<-results[grepneun,]
neun$Category<-unlist(lapply(strsplit(neun$Category, split="L2_0"), '[[', 1))

neunenrich<-neunpval<-neunherit<-neunse<-matrix(0, nrow=length(neun$Category), ncol=length(conditions))
rownames(neunenrich)<-rownames(neunpval)<-rownames(neunherit)<-rownames(neunse)<-neun$Category
colnames(neunenrich)<-colnames(neunpval)<-colnames(neunherit)<-colnames(neunse)<-conditions

for(i in 1:length(conditions)){
    condition<-conditions[i]
    results<-read.table(paste0(condition,"_results.results"), header=T)
    #grepneun<-grep("NeuN", results$Category)
    neun<-results[grepneun,]
    neun$Category<-unlist(lapply(strsplit(neun$Category, split="L2_0"), '[[', 1))
    neunenrich[,i]<-neun$Enrichment
    neunpval[,i]<-neun$Enrichment_p
    neunherit[,i]<-neun$Prop._h2
    neunse[,i]<-neun$Enrichment_std_error
    print(condition)
}

save(neunpval, neunenrich, neunherit, neunse, conditions,
     file="cell_type_Glu_GABA_enrichment_LDSC.rda")

load("cell_type_Glu_GABA_enrichment_LDSC.rda")

neunpvalmelt<-melt(neunpval)
neunenrichmelt<-melt(neunenrich)
neunsemelt<-melt(neunse)
colnames(neunpvalmelt)<-colnames(neunenrichmelt)<-colnames(neunsemelt)<-c('celltype','conditions','value')
neunpvalmelt$FDR<-p.adjust(neunpvalmelt$value, "BH")

cellname<-c('SPARK_EUR_iPSYCH_PGC_Meta','adhd2','mdd.2019','scz3','bp2019')

pdf(file=paste0("cell_type_Glu_GABA_enrichment_disorder.pdf"), width=8, height=5)
barpl<-ggplot(neunpvalmelt, aes(x=factor(conditions,levels=cellname),y=-log10(value),fill=celltype))
barpl + geom_bar(stat="identity",position="dodge") + theme_bw() + ggtitle("P-value Enrichment") + scale_fill_brewer(palette="Paired") +
  theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(), axis.text.x=element_text(angle=90,hjust=1)) + ylab("-log10(P-value)") + 
  geom_hline(yintercept=-log10(0.05),lty=2,color="red")

barp2<-ggplot(neunpvalmelt, aes(x=factor(conditions,levels=cellname),y=-log10(FDR),fill=celltype))
barp2 + geom_bar(stat="identity",position="dodge") + theme_bw() + ggtitle("FDR Enrichment") + scale_fill_brewer(palette="Paired") +
  theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(), axis.text.x=element_text(angle=90,hjust=1)) + ylab("-log10(FDR)") + 
  geom_hline(yintercept=-log10(0.05),lty=2,color="red") 

barp3<-ggplot(neunenrichmelt, aes(x=factor(conditions,levels=cellname),y=value,fill=celltype))
barp3 + geom_bar(stat="identity",position="dodge",colour="white") + theme_bw() + ggtitle("Enrichment_score") + ylab("Heritability Enrichment") + scale_fill_brewer(palette="Paired") +
  theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(), axis.text.x=element_text(angle=90,hjust=1)) + geom_hline(yintercept=1,lty=2,color="red") + 
  geom_errorbar(aes(ymin=value-neunsemelt$value, ymax=value+neunsemelt$value), width=.2,position=position_dodge(.9))

dev.off()
