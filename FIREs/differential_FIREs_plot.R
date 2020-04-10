#!/usr/bin/env Rscript
library(ggplot2)
library(reshape2)
library(gProfileR)
library("optparse")
 options(stringsAsFactors = FALSE)
option_list<-list(
    make_option(c("-f", "--fname"), type="character", default=NULL, 
            help="dataset  file name", metavar="character"),
    make_option(c("-o", "--oname"), type="character", default=NULL, 
            help="output file name", metavar="character")
)

opt_parser<-OptionParser(option_list=option_list)
opt<-parse_args(opt_parser)

#which produces a list opt that contains all the arguments sorted by order of appearance in option_list
# and which can be called by their names as declared in this object: opt$file and opt$out. 
#Then, managing null arguments is performed as follows:
if (is.null(opt$fname)){
    print_help(opt_parser)
    stop("At least one argument must be supplied (fire file).n", call.=FALSE)
}
#in which the function print_help print the help page of the option list as declared in the object option_list.
## program...
fire<-read.table(opt$fname, header=TRUE)
fire<-na.omit(fire)
fire$colname<-c("gray")
fire$pchname<-c(1)
fire$colname[fire$NeuNpos_FIRES>qnorm(0.975) & fire$NeuNneg_FIRES<qnorm(0.9)]<-"red"
fire$pchname[fire$NeuNpos_FIRES>qnorm(0.975) & fire$NeuNneg_FIRES<qnorm(0.9)]<-16
fire$colname[fire$NeuNpos_FIRES<qnorm(0.9) & fire$NeuNneg_FIRES>qnorm(0.975)]<-"blue"
fire$pchname[fire$NeuNpos_FIRES<qnorm(0.9) & fire$NeuNneg_FIRES>qnorm(0.975)]<-18
pdf(paste0(opt$oname,"_diffFIREs_plot.pdf"))
p<-plot(fire$NeuNpos_FIRES,fire$NeuNneg_FIRES,col=fire$colname,pch=fire$pchname,xlab="NeuNpos FIREs",ylab="NeuNneg FIRES")
legend("topleft",pch=c(16,18,1),col=c("red","blue","gray"),c("neuornal FIREs","glial FIREs","common FIREs"))
abline(h=qnorm(0.975),lty=2)
abline(h=qnorm(0.9),lty=2)
abline(v=qnorm(0.975),lty=2)
abline(v=qnorm(0.9),lty=2)
print(p)
dev.off()
pdf(paste0(opt$oname,"_diffFIREs_plot2.pdf"))
ggplot(fire, aes(x=NeuNpos_FIRES, y=NeuNneg_FIRES, shape=colname, color=colname)) +
    geom_point()+ scale_color_manual(values=c("blue", "gray", "red"))+ 
    geom_hline(yintercept=qnorm(0.975), linetype="dashed")+
    geom_hline(yintercept=qnorm(0.9), linetype="dashed")+
    geom_vline(xintercept = qnorm(0.975), linetype="dashed")+
    geom_vline(xintercept = qnorm(0.9), linetype="dashed")
dev.off()
#If the entire script is saved in a file called yasrs.R

#Rscript --vanilla yasrs.R  -f iris.txt -o out.txt