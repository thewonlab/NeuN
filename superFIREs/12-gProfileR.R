setwd("/nas/longleaf/home/hubenxia/project/FIRE/neuronvsglia/superFIRE")
library(gProfileR)
options(stringsAsFactors = FALSE)
posgene<-read.table("NeuNposgene_superFIRE_H3K27ac.txt",header=T)[,1] ##superFIREs overlapping with gene promoters
neggene<-read.table("NeuNneggene_superFIRE_H3K27ac.txt",header=T)[,1]

GO<-function(gene_list,flag){
    goresult <- gprofiler(gene_list,organism = "hsapiens",ordered_query=F,significant=T,
            max_p_value=0.1,min_set_size=15,max_set_size=600,
            min_isect_size=5,correction_method="gSCS",
            hier_filtering="strong",include_graph=T,src_filter="GO")
    write.table(file=paste0("goresult_",flag,"_superFIRE_H3K27ac.txt"),goresult,quote=F,row.names=F,sep="\t")
    return(goresult)
}
pos_GO<-GO(posgene,"NeuNposgene")
neg_GO<-GO(neggene,"NeuNneggene")
library(ggplot2)
GO<-function(input,filename){
    GO_order<-input[order(input$p.value),][c(1:10),]
    pdf(paste0("GO_terms_for_",filename,".pdf"))
    p<-ggplot(data=GO_order, aes(y=-log10(p.value), x=reorder(term.name,-log10(p.value)))) +
    geom_bar(stat="identity", width=0.5)+labs(title=paste0("GO terms for ",filename),y="-log10(p value)", x = "Term names")+ coord_flip()+
    theme_classic()+ geom_hline(yintercept =-log10(0.05),color = "red")
    print(p)
    dev.off()
}

GO(neg_GO,'NeuNneg_superFIRE_H3K27ac')
GO(pos_GO,'NeuNpos_superFIRE_H3K27ac')
