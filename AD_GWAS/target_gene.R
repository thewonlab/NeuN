#!/usr/bin/env Rscript
library(optparse)
options(stringsAsFactors = FALSE)
library(biomaRt)
library(reshape)
library(WGCNA)
library(ggplot2)
library(gridExtra)
library(GenomicRanges)

option_list<-list(
    make_option(c("-a", "--anno"), type="character", default=NULL, 
            help="anno file name", metavar="character"),
    make_option(c("-i", "--input"), type="character", default=NULL, 
            help="input file name", metavar="character"),
    make_option(c("-d", "--disease"), type="character", default=NULL, 
            help="disease file name", metavar="character"),
    make_option(c("-o", "--outdir"), type="character", default=NULL, 
            help="outdir file name", metavar="character")
)

opt_parser<-OptionParser(option_list=option_list)
opt<-parse_args(opt_parser)

if (is.null(opt$anno)){
    print_help(opt_parser)
    stop("At least one argument must be supplied (anno file).n", call.=FALSE)
}
if (is.null(opt$input)){
    print_help(opt_parser)
    stop("At least one argument must be supplied (input file).n", call.=FALSE)
}
if (is.null(opt$disease)){
    print_help(opt_parser)
    stop("At least one argument must be supplied (disease file).n", call.=FALSE)
}

if (is.null(opt$outdir)){
    print_help(opt_parser)
    stop("At least one argument must be supplied (outdir file).n", call.=FALSE)
}

setwd(opt$outdir)
load(opt$anno) # Annotation for all genes 
pcgene = unique(geneAnno1[geneAnno1$gene_biotype=="protein_coding","ensembl_gene_id"])  # protein-coding genes
diseasename = c(opt$disease)

celltype = c("glia","neuron")
fdrthr= 0.05
mhcrange = GRanges(seqnames=6, IRanges(25000000,35000000)) # MHC region
getout<-function(filename,celltype,mhcrange,diseasename){
    diseasemat = read.table(paste0(filename,"_",celltype,".genes.out"), header=T) # reading MAGMA output files
    diseasepc = diseasemat[diseasemat$CHR!="X" & diseasemat$GENE %in% pcgene, ] # only get protein-coding genes
    # selecting out genes that are located in MHC region
    diseaseranges = GRanges(seqnames=diseasepc$CHR, IRanges(diseasepc$START, diseasepc$STOP), gene=diseasepc$GENE)
    olap = findOverlaps(diseaseranges, mhcrange)
    mhcgene = diseaseranges[queryHits(olap)]$gene
    diseasepc = diseasepc[!(diseasepc$GENE %in% mhcgene),] # remove MHC genes
    diseasepc$fdr = p.adjust(diseasepc$P, "BH") # calculate FDR from P-value
    selectgenes = diseasepc[diseasepc$fdr<fdrthr, "GENE"] # select ouf significantly associated genes
    diseasefdr = selectgenes
    genename = geneAnno1[geneAnno1$ensembl_gene_id %in% diseasefdr, "hgnc_symbol"]
    genename = genename[genename!=""]

    backname = geneAnno1[geneAnno1$ensembl_gene_id %in% diseasepc$GENE, "hgnc_symbol"]
    backname = backname[backname!=""]
    print(paste(filename,length(diseasefdr)))
    save(diseasefdr, genename,diseasepc, file=paste0(opt$outdir,diseasename,"_",celltype,"_fdrthrgenes_adjusted.rda"))
    write.table(diseasefdr, file=paste0(opt$outdir,diseasename,"_",celltype,'_geneid.txt'),quote=F,row.names=F,col.names=F)
    write.table(genename, file=paste0(opt$outdir,diseasename,"_",celltype,'_genename.txt'),quote=F,row.names=F,col.names=F)
    write.table(diseasepc$GENE, file=paste0(opt$outdir,diseasename,"_",celltype,'_backgeneid.txt'),quote=F,row.names=F,col.names=F)
    write.table(backname, file=paste0(opt$outdir,diseasename,"_",celltype,'_backgenename.txt'),quote=F,row.names=F,col.names=F)
}

getout(opt$input,celltype[1],mhcrange,diseasename)
getout(opt$input,celltype[2],mhcrange,diseasename)