#!/usr/bin/env Rscript
library(ggplot2)
library(reshape2)
library(gProfileR)
library("optparse")
library(GenomicRanges)
library(data.table)
library(stringr)
library(stats4)
library(BiocGenerics)
library(parallel)
library(S4Vectors)
library(IRanges)
library(GenomeInfoDb)
options(stringsAsFactors = FALSE)
option_list<-list(
    make_option(c("-f", "--fname"), type="character", default=NULL, 
            help="dataset  file name", metavar="character"),
    make_option(c("-a", "--anno"), type="character", default=NULL, 
            help="anno  file name", metavar="character"),
    make_option(c("-o", "--oname"), type="character", default=NULL, 
            help="pdf file name", metavar="character")
)

opt_parser<-OptionParser(option_list=option_list)
opt<-parse_args(opt_parser)

if (is.null(opt$fname)){
    print_help(opt_parser)
    stop("At least one argument must be supplied (input file).n", call.=FALSE)
}
if (is.null(opt$anno)){
    print_help(opt_parser)
    stop("At least one argument must be supplied (anno file).n", call.=FALSE)
}
load(opt$anno)
geneAnno1<-geneAnno1[geneAnno1$gene_biotype=="protein_coding" & geneAnno1$hgnc_symbol!="",]
diseasemat<-read.table(opt$fname, header=T)
diseasemat <- diseasemat[diseasemat$GENE %in% geneAnno1$ensembl_gene_id, ]
diseaseranges <- GRanges(seqnames=diseasemat$CHR, IRanges(diseasemat$START, diseasemat$STOP), ensg=diseasemat$GENE, P=diseasemat$P)
MHCrange <- GRanges(seqnames=6, IRanges(25000000,35000000)) 
olap <- findOverlaps(diseaseranges, MHCrange)
diseasemhc <- diseaseranges[queryHits(olap)]
mhcgene <- diseasemhc$ensg
diseasemat <- diseasemat[!(diseasemat$GENE %in% mhcgene), ]
backgroundensg <- diseasemat$GENE
backgroundensg <- geneAnno1[geneAnno1$ensembl_gene_id %in% backgroundensg, "hgnc_symbol"]

diseasemat<-merge(diseasemat,geneAnno1,by.x='GENE',by.y='ensembl_gene_id')
queryensg <- diseasemat[order(diseasemat$P),"hgnc_symbol"]

write.table(file=paste0("bg_",opt$oname,".txt"),backgroundensg,quote=F,col.names=F,row.names=F,sep="\t")
write.table(file=paste0("qg_",opt$oname,".txt"),queryensg,quote=F,col.names=F,row.names=F,sep="\t")

GO<-function(input,filename,bg){
    goresult<-gprofiler(input,  organism="hsapiens",  ordered_query=T, significant=T,  max_p_value=0.1,
    min_set_size=15,max_set_size=600, min_isect_size=5,custom_bg=bg, correction_method="gSCS",  hier_filtering="moderate", 
     include_graph=T, src_filter="GO")
    write.table(file=paste0("goresult_",filename,".txt"),goresult,quote=F,row.names=F,sep="\t")
    return(goresult)
}

result<-GO(queryensg,opt$oname,backgroundensg)

GO<-function(input,filename){
    GO_order<-input[order(input$p.value),][c(1:10),]
    pdf(paste0("GO_terms_for_",filename,".pdf"))
    p<-ggplot(data=GO_order, aes(y=-log10(p.value), x=reorder(term.name,-log10(p.value))))+ geom_hline(yintercept =-log10(0.05),color = "red") +
    geom_bar(stat="identity", width=0.5)+labs(title=paste0("GO terms for ",filename),y="-log10(p value)", x = "Term names")+ coord_flip()+
    theme_classic()
    print(p)
    dev.off()
}

GO(result,opt$oname)