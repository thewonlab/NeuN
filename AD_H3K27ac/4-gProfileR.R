library(gProfileR)
pos_hypo<-read.table("pos_hypo.genelist.txt",header=T)
neg_hypo<-read.table("neg_hypo.genelist.txt",header=T)
pos_hyper<- read.table("pos_hyper.genelist.txt", header=T,sep="\t")
neg_hyper<-read.table("neg_hyper.genelist.txt", header=T,sep="\t")

go<-function(gene,sample_name){
    gene_list<-as.vector(gene[,1])
    goresult = gprofiler(gene_list,organism = "hsapiens",ordered_query=F,significant=T,
            correction_method="gSCS",
            hier_filtering="strong",include_graph=T,src_filter="GO")
    write.table(file=paste0("goresult_",sample_name,"_Alzheimer.txt"),goresult,quote=F,row.names=F,sep="\t")
}
go(pos_hypo,"pos_hypo")
go(neg_hypo,"neg_hypo")
go(pos_hyper,"pos_hyper")
go(neg_hyper,"neg_hyper")
