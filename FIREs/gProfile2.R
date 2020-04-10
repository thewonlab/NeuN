#!/usr/bin/env Rscript
library(ggplot2)
library(gprofiler2)
library(forcats)
library(optparse)
options(stringsAsFactors=F)
option_list<-list(
    make_option(c("-f", "--fname"), type="character", default=NULL, 
            help="dataset  file name", metavar="character"),
    make_option(c("-o", "--oname"), type="character", default=NULL, 
            help="pdf file name", metavar="character")
)

opt_parser<-OptionParser(option_list=option_list)
opt<-parse_args(opt_parser)

if (is.null(opt$fname)){
    print_help(opt_parser)
    stop("At least one argument must be supplied (input file).n", call.=FALSE)
}
if (is.null(opt$oname)){
    print_help(opt_parser)
    stop("At least one argument must be supplied (oname file).n", call.=FALSE)
}
gene<-read.table(opt$fname,header=F,sep="\t")[,1]
print(head(gene))
GO<-function(input,filename){
    goresult <- gost(input,organism = "hsapiens",ordered_query=F,significant=T,evcodes = TRUE,
            user_threshold=0.05, correction_method="fdr",sources=c("GO"))$result
    goresult <- data.frame(term_id=goresult$term_id,p_value=goresult$p_value,term_name=goresult$term_name,gene=goresult$intersection,size=goresult$term_size)
    goresult <- goresult[goresult$size>15 & goresult$size<600,]
    write.table(file=paste0("goresult_",filename,".txt"),goresult,quote=F,row.names=F,sep="\t")
    write.table(file=paste0("goresult_",filename,"_REViGO.txt"),goresult[,1:2],quote=F,row.names=F,sep="\t")
    return(goresult)
}

result<-GO(gene,opt$oname)

gplot<-function(input,filename){
    GO_order<-input[order(input$p_value),][c(1:10),]
    pdf(paste0("GO_terms_for_",filename,".pdf"))
    p<-ggplot(data=GO_order, aes(y=-log10(p_value), x=reorder(term_name,-log10(p_value)))) +
    geom_bar(stat="identity", width=0.5)+labs(title=paste0("GO terms for ",filename),y="-log10(p value)", x = "Term names")+ coord_flip()+
    theme_classic()+ geom_hline(yintercept =-log10(0.05),color = "red")
    print(p)
    dev.off()
}
gplot(result,opt$oname)