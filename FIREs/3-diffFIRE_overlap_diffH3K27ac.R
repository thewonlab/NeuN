#!/usr/bin/env Rscript
library(ggplot2)
library(reshape2)
library("optparse")
options(stringsAsFactors=F)
library(stats4)
library(BiocGenerics)
library(parallel)
library(S4Vectors)
library(IRanges)
library(GenomeInfoDb)
library(GenomicRanges)
library(dplyr)

option_list<-list(
    make_option(c("-f", "--FIRE"), type="character", default=NULL, 
            help="FIRE  file name", metavar="character"),
    make_option(c("-n", "--neuronalpeak"), type="character", default=NULL, 
            help="neuronal H3K27ac peak file name", metavar="character"),
    make_option(c("-g", "--glialpeak"), type="character", default=NULL, 
            help="glial H3K27ac peak file name", metavar="character"),
    make_option(c("-N", "--NeuronOut"), type="character", default=NULL, 
            help="neuronal output file name", metavar="character"),
    make_option(c("-G", "--GliaOut"), type="character", default=NULL, 
            help="glial output file name", metavar="character")
)

opt_parser<-OptionParser(option_list=option_list)
opt<-parse_args(opt_parser)

if (is.null(opt$FIRE)){
    print_help(opt_parser)
    stop("At least one argument must be supplied (FIRE file).n", call.=FALSE)
}
if (is.null(opt$glialpeak)){
    print_help(opt_parser)
    stop("At least one argument must be supplied (glial peak file).n", call.=FALSE)
}
if (is.null(opt$neuronalpeak)){
    print_help(opt_parser)
    stop("At least one argument must be supplied (neuronal peak file).n", call.=FALSE)
}

#### load data
fire<-read.table(opt$FIRE, header=T)
neuron<-read.table(opt$neuronalpeak,sep="\t",header=T,quote="")
glia<-read.table(opt$glialpeak,sep="\t",header=T,quote="")
#### select differential FIREs
firediff.NeuNpos<-fire[fire$NeuNpos_FIRES>qnorm(0.975) & fire$NeuNneg_FIRES<qnorm(0.9),]
firediff.NeuNneg<-fire[fire$NeuNpos_FIRES<qnorm(0.9) & fire$NeuNneg_FIRES>qnorm(0.975),]

pos_fire<-nrow(firediff.NeuNpos)
neg_fire<-nrow(firediff.NeuNneg)

firanges.NeuNpos<-GRanges(seqnames=firediff.NeuNpos$chr, IRanges(firediff.NeuNpos$start, firediff.NeuNpos$end))
firanges.NeuNneg<-GRanges(seqnames=firediff.NeuNneg$chr, IRanges(firediff.NeuNneg$start, firediff.NeuNneg$end))

#####################
overalapping<-function(input,target){
    olap<-findOverlaps(input, target)
    result<-data.frame(input[queryHits(olap)])[,c(1:3)]
    return(nrow(unique(result)))
}
neuronranges<-GRanges(seqnames=neuron$Chrom, IRanges(neuron$Start, neuron$End))
gliaranges<-GRanges(seqnames=glia$Chrom, IRanges(glia$Start, glia$End))

output<-function(total,neuron,glia){
    fire_neuron<-overalapping(total,neuron)
    fire_glia<-overalapping(total,glia)
    total_num<-nrow(data.frame(total))
    print(total)
    print(total_num)
    fire_neuron_other<-total_num-fire_neuron
    fire_glia_other<-total_num-fire_glia
    peak<-c(fire_neuron,fire_glia)
    other<-c(fire_neuron_other,fire_glia_other)
    name<-c('neuron','glia')
    fire_peak<-data.frame(type=name,peak=peak,other=other)
    result<-melt(fire_peak)
    return(result)
}


#######################NeuNpos overlapping with differential H3K27ac peaks
result_pos<-output(firanges.NeuNpos,neuronranges,gliaranges)
print(result_pos)
pdf(paste0('geom_bar_pos_diffFIRE_',opt$NeuronOut,'_H3K27ac_peak.pdf'), width=5, height=5)
ggplot(result_pos,aes(x=type,y=value,fill=variable))+
geom_bar(stat="identity",position ="fill")+
labs(title="proportion of diffFIRE overlapping with peak",y="proportion", x = "Hi-C samples")+
theme_bw() + theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank())
dev.off()

################# NeuNneg overlapping with differential H3K27ac peaks
result_neg<-output(firanges.NeuNneg,neuronranges,gliaranges)
print(result_neg)
pdf(paste0('geom_bar_neg_diffFIRE_',opt$GliaOut,'_H3K27ac_peak.pdf'), width=5, height=5)
ggplot(result_neg,aes(x=type,y=value,fill=variable))+
geom_bar(stat="identity",position ="fill")+
labs(title="proportion of diffFIRE overlapping with peak",y="proportion", x = "Hi-C samples")+
theme_bw() + theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank())
dev.off()
