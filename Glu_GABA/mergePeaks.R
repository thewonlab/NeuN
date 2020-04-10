#!/usr/bin/env Rscript
library(optparse)
options(stringsAsFactors=F)
library(GenomicRanges)
library(stats4)
library(BiocGenerics)
library(parallel)
library(S4Vectors)
library(IRanges)
library(GenomeInfoDb)
library(rtracklayer)
library(GenomicRanges)

#merge peaks from multiple replicates

option_list<-list(
    make_option(c("-d", "--dir"), type="character", default=NULL, 
            help="directory name", metavar="character"),
    make_option(c("-p", "--peak"), type="character", default=NULL, 
            help="peak file name", metavar="character"),
    make_option(c("-o", "--output"), type="character", default=NULL, 
            help="merged file name", metavar="character")
)

opt_parser<-OptionParser(option_list=option_list)
opt<-parse_args(opt_parser)

if (is.null(opt$dir)){
    print_help(opt_parser)
    stop("At least one argument must be supplied (directory file).n", call.=FALSE)
}
if (is.null(opt$peak)){
    print_help(opt_parser)
    stop("At least one argument must be supplied (peak file).n", call.=FALSE)
}
if (is.null(opt$output)){
    print_help(opt_parser)
    stop("At least one argument must be supplied (output file).n", call.=FALSE)
}

setwd(opt$dir)
if(opt$peak == 'GLU'){
    peak_files <- list.files(pattern = "*.GLU_peaks.broadPeak")
}else{
    peak_files <- list.files(pattern = "*.SOX_peaks.broadPeak")
}
peak_files
peak_grangeslist <-list()
for (a in peak_files){
    input<-read.table(a,sep="\t",header=F)
    result<-GRanges(seqnames=input[,1],IRanges(input[,2],input[,3]))
    mcols(result)<-input[,c(4:5,7:9)]
    peak_grangeslist <-append(peak_grangeslist,list(result))
}

peak_grangeslist <- GRangesList(peak_grangeslist)
peak_grangeslist
peak_coverage <- coverage(peak_grangeslist)
peak_coverage
covered_ranges <- slice(peak_coverage, lower=4, rangesOnly=T)
covered_ranges
covered_granges <- GRanges(covered_ranges)
covered_granges
result<-reduce(covered_granges, min.gapwidth=1)
export(result, opt$output)