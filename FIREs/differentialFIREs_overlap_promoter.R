#!/usr/bin/env Rscript
library("optparse")
options(stringsAsFactors = FALSE)
library(IRanges)
library(GenomeInfoDb)
library(GenomicRanges)
option_list<-list(
    make_option(c("-p", "--pname"), type="character", default=NULL, 
            help="promoter  file name", metavar="character"),
    make_option(c("-f", "--fname"), type="character", default=NULL, 
            help="FIREs file name", metavar="character"),
    make_option(c("-g", "--gname"), type="character", default=NULL, 
            help="annotation file name", metavar="character"),
    make_option(c("-o", "--oname"), type="character", default=NULL, 
            help="output file name", metavar="character")
)

opt_parser<-OptionParser(option_list=option_list)
opt<-parse_args(opt_parser)

if (is.null(opt$pname)){
    print_help(opt_parser)
    stop("At least one argument must be supplied (pname file).n", call.=FALSE)
}
if (is.null(opt$fname)){
    print_help(opt_parser)
    stop("At least one argument must be supplied (fname file).n", call.=FALSE)
}
if (is.null(opt$gname)){
    print_help(opt_parser)
    stop("At least one argument must be supplied (gname file).n", call.=FALSE)
}
if (is.null(opt$oname)){
    print_help(opt_parser)
    stop("At least one argument must be supplied (oname file).n", call.=FALSE)
}
fire<-read.table(opt$fname, header=T,sep="\t")
promoter<-read.table(opt$pname)
promoteranges<-GRanges(seqnames=promoter[,1], IRanges(promoter[,2],promoter[,3]), gene=promoter[,6])
load(opt$gname)
firanges<-GRanges(seqnames=fire$chr, IRanges(fire$start,fire$end))
overalapping<-function(input,target,anno,flag){
    olap <- findOverlaps(input, target)
    #result <- input[queryHits(olap)]
    #mcols(result) <- mcols(target[subjectHits(olap)])
    gene <- unique(mcols(target[subjectHits(olap)]))[,1]
    write.table(file=paste0(flag,'_diffFIRE_promoter_geneID.txt'),gene,sep="\t",quote=F,col.names=F,row.names=F)
    gene <- anno[(anno[,1] %in% gene) & anno[,2]!='',2]
    write.table(file=paste0(flag,'_diffFIRE_promoter_gene.txt'),gene,sep="\t",quote=F,col.names=F,row.names=F)
}

overalapping(firanges,promoteranges,geneAnno1,opt$oname)