#!/usr/bin/env Rscript
library(ggplot2)
library(reshape2)
library(gProfileR)
library("optparse")
 options(stringsAsFactors = FALSE)
option_list<-list(
    make_option(c("-f", "--fname"), type="character", default=NULL, 
            help="dataset  file name", metavar="character"),
    make_option(c("-N", "--Nname"), type="character", default=NULL, 
            help="neuornal FIREs file name", metavar="character"),
    make_option(c("-G", "--Gname"), type="character", default=NULL, 
            help="glial FIREs file name", metavar="character")
)

opt_parser<-OptionParser(option_list=option_list)
opt<-parse_args(opt_parser)

#which produces a list opt that contains all the arguments sorted by order of appearance in option_list
# and which can be called by their names as declared in this object: opt$file and opt$out. 
#Then, managing null arguments is performed as follows:
if (is.null(opt$fname)){
    print_help(opt_parser)
    stop("At least one argument must be supplied (input file).n", call.=FALSE)
}
#in which the function print_help print the help page of the option list as declared in the object option_list.
## program...
fire<-read.table(opt$fname, header=TRUE)
firediff.NeuNpos<-fire[fire$NeuNpos_FIRES>qnorm(0.975) & fire$NeuNneg_FIRES<qnorm(0.9),]
NeuNpos<-fire[fire$NeuNpos_FIRES>qnorm(0.975) & fire$NeuNneg_FIRES<qnorm(0.9),c(1:3)]
nrow(firediff.NeuNpos)
firediff.NeuNneg<-fire[fire$NeuNpos_FIRES<qnorm(0.9) & fire$NeuNneg_FIRES>qnorm(0.975),]
NeuNneg<-fire[fire$NeuNpos_FIRES<qnorm(0.9) & fire$NeuNneg_FIRES>qnorm(0.975),c(1:3)]
nrow(firediff.NeuNneg)

write.table(firediff.NeuNpos,file=paste0(opt$Nname,"_diffFIREs.bed"),sep="\t",quote=F,row.names=F)
write.table(firediff.NeuNneg,file=paste0(opt$Gname,"_diffFIREs.bed"),sep="\t",quote=F,row.names=F)
write.table(NeuNpos,file=paste0(opt$Nname,"_diffFIREs_GREAT.bed"),sep="\t",quote=F,row.names=F,col.names=F)
write.table(NeuNneg,file=paste0(opt$Gname,"_diffFIREs_GREAT.bed"),sep="\t",quote=F,row.names=F,col.names=F)
#If the entire script is saved in a file called yasrs.R

#Rscript --vanilla yasrs.R  -f iris.txt -o out.txt
