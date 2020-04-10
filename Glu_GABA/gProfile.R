#!/usr/bin/env Rscript
library(ggplot2)
library(gProfileR)
library(optparse)
options(stringsAsFactors=F)
option_list<-list(
    make_option(c("-f", "--fname"), type="character", default=NULL, 
            help="fname file name", metavar="character"),
    make_option(c("-o", "--outname"), type="character", default=NULL, 
            help="outname file name", metavar="character")
)

opt_parser<-OptionParser(option_list=option_list)
opt<-parse_args(opt_parser)
if (is.null(opt$fname)){
    print_help(opt_parser)
    stop("At least one argument must be supplied (fname file).n", call.=FALSE)
}
if (is.null(opt$outname)){
    print_help(opt_parser)
    stop("At least one argument must be supplied (outname file).n", call.=FALSE)
}
gene<-read.table(opt$fname,header=F)[,1]

GO<-function(input,filename){
    goresult <- gprofiler(input,organism = "hsapiens",ordered_query=F,significant=T,
            max_p_value=0.1,min_set_size=5, max_set_size=1000,
            min_isect_size=5,correction_method="gSCS",
            hier_filtering="strong",include_graph=T,src_filter="GO")
    write.table(file=paste0("goresult_",filename,"_loops.txt"),goresult,quote=F,row.names=F,sep="\t")
    return(goresult)
}

result<-GO(gene,opt$outname)

gplot<-function(input,filename){
    GO_order<-input[order(input$p.value),][c(1:10),]
    pdf(paste0("GO_terms_for_",filename,"_Loops.pdf"))
    p<-ggplot(data=GO_order, aes(y=-log10(p.value), x=reorder(term.name,-log10(p.value)))) +
    geom_bar(stat="identity", width=0.5)+labs(title=paste0("GO terms for ",filename),y="-log10(p value)", x = "Term names")+ coord_flip()+
    theme_classic()
    print(p)
    dev.off()
}
gplot(result,opt$outname)