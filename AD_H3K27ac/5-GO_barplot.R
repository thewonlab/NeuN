library(ggplot2)
pos_hypo<-read.table("goresult_pos_hypo_Alzheimer.txt",sep="\t",header=T)
neg_hyper<-read.table("goresult_neg_hyper_Alzheimer.txt",sep="\t",header=T)

bar<-function(GO,sample_name){
    GO_order<-GO[order(GO$p.value),]
    pdf(paste0("GO_terms_for_",sample_name,".pdf"))
    p<-ggplot(data=GO_order, aes(y=-log10(p.value), x=reorder(term.name,-log10(p.value)))) +
    geom_bar(stat="identity", width=0.5)+labs(title=paste0("GO_terms_for_",sample_name),y="-log10(p value)", x = "Term names")+ coord_flip()+theme_classic()
    print(p)  ####must add this code
    dev.off()
}

bar(pos_hypo,"pos_hypo")
bar(neg_hyper,"neg_hyper")