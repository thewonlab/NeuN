#!/usr/bin/env Rscript
library(optparse)
options(stringsAsFactors = FALSE)
library(biomaRt)
library(reshape)
library(WGCNA)
library(ggplot2)
library(gridExtra)
library(GenomicRanges)
library(stats4)
library(BiocGenerics)
library(parallel)
library(S4Vectors)
library(IRanges)
library(GenomeInfoDb)
library(reshape2)
library(dplyr)
library(corrplot)

option_list<-list(
    make_option(c("-S", "--GABA"), type="character", default=NULL, 
            help="GABA file name", metavar="character"),
    make_option(c("-B", "--Glu"), type="character", default=NULL, 
            help="Glu file name", metavar="character"),
    make_option(c("-s", "--scRNAseq"), type="character", default=NULL, 
            help="scRNAseq file name", metavar="character"),
    make_option(c("-N", "--Neuron"), type="character", default=NULL, 
            help="Neuron name", metavar="character"),
    make_option(c("-G", "--Glia"), type="character", default=NULL, 
            help="Glia name", metavar="character"),
    make_option(c("-d", "--disease"), type="character", default=NULL, 
            help="disease name", metavar="character"),
    make_option(c("-o", "--outdir"), type="character", default=NULL, 
            help="outdir file name", metavar="character")
)

opt_parser<-OptionParser(option_list=option_list)
opt<-parse_args(opt_parser)

if (is.null(opt$GABA)){
    print_help(opt_parser)
    stop("At least one argument must be supplied (GABA file).n", call.=FALSE)
}
if (is.null(opt$Glu)){
    print_help(opt_parser)
    stop("At least one argument must be supplied (Glu file).n", call.=FALSE)
}
if (is.null(opt$scRNAseq)){
    print_help(opt_parser)
    stop("At least one argument must be supplied (scRNAseq file).n", call.=FALSE)
}
if (is.null(opt$Neuron)){
    print_help(opt_parser)
    stop("At least one argument must be supplied (Neuron name).n", call.=FALSE)
}
if (is.null(opt$Glia)){
    print_help(opt_parser)
    stop("At least one argument must be supplied (Glia name).n", call.=FALSE)
}
if (is.null(opt$disease)){
    print_help(opt_parser)
    stop("At least one argument must be supplied (disease name).n", call.=FALSE)
}
if (is.null(opt$outdir)){
    print_help(opt_parser)
    stop("At least one argument must be supplied (outdir file).n", call.=FALSE)
}

load(opt$GABA)
neuron = genename
load(opt$Glu)
glia = genename

cellexp = read.table(opt$scRNAseq, header=T)
rownames(cellexp)<-cellexp[,1]
cellexp<-cellexp[,-1]
datExpr<-scale(cellexp,center=T, scale=F) ###
theme_set(theme_bw())

gplot<-function(datExpr,targetgene,disease,filename){
    exprdat = apply(datExpr[match(targetgene, 
        rownames(datExpr)),],2,mean,na.rm=T) 
    dat = data.frame(Group=disease, cell=names(exprdat), Expr=exprdat) 
    dat$cell = gsub("fetal.quiescent","fetal-quiescent",dat$cell)
    dat$cell = gsub("fetal.replicating","fetal-replicating",dat$cell)
    dat$celltype = unlist(lapply(strsplit(dat$cell, split="[.]"),'[[',1))
    dat = (dat[(dat$celltype %in% c("Ex1","Ex2","Ex3","Ex4","Ex5","Ex6","Ex7","Ex8","In1","In2","In3","In4","In5","In6","In7","In8")),])
    dat$celltype = factor(dat$celltype, levels=c("Ex1","Ex2","Ex3","Ex4","Ex5","Ex6","Ex7","Ex8","In1","In2","In3","In4","In5","In6","In7","In8"))

    #pdf(file=paste0("singlecell_expression_",filename,".pdf"))
    #p<-ggplot(dat,aes(x=celltype, y=Expr, fill=celltype)) +  
    #   ylab("Normalized expression") + xlab("") + geom_violin() + 
    #   theme(axis.text.x=element_text(angle = 90, hjust=1)) + theme_classic()+ 
    #   theme(legend.position="none") +  ggtitle(paste0("Cellular expression profiles of ",filename)) 
    #print(p)
    #dev.off() 
}
print(paste0(opt$Neuron,opt$disease,'gene'))
gplot(datExpr,neuron,opt$disease,paste0(opt$Neuron,opt$disease,'gene'))
gplot(datExpr,glia,opt$disease,paste0(opt$Glia,opt$disease,'gene'))

summarySE<-function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
    library(plyr)
  # New version of length which can handle NA's: if na.rm==T, don't count them
    length2<-function(x, na.rm=FALSE) {
        if (na.rm){
            sum(!is.na(x))
        }else{
            length(x)
        }
    }
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
    datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
    #Rename the "mean" column    
    datac <- rename(datac, c("mean" = measurevar))
    datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
    #Confidence interval multiplier for standard error
    #Calculate t-statistic for confidence interval: 
    #e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
    ciMult <- qt(conf.interval/2 + .5, datac$N-1)
    datac$ci <- datac$se * ciMult
    return(datac)
}


target<-function(data1,data2,name1,name2){
    targetname<-c(name1, name2)
    targetgene<-vector(length=2, mode="list")
    targetgene <- list(data1,data2)
    

    exprdat<-vector(mode="list", length=length(targetgene))

    for(i in 1:length(targetgene)){
        exprdat[[i]]<-apply(datExpr[match(targetgene[[i]], rownames(datExpr)),],2,mean,na.rm=T)
    }

    dat<-c()
    for(i in 1:length(targetgene)){
        datframe = data.frame(Group=targetname[i], cell=names(exprdat[[i]]), Expr=exprdat[[i]])
        dat = rbind(dat, datframe)
    }

    dat$cell<-gsub("Fetal.quiescent","Fetal-quiescent",dat$cell)
    dat$cell<-gsub("Fetal.replicating","Fetal-replicating",dat$cell)
    dat$celltype<- unlist(lapply(strsplit(dat$cell, split="[.]"),'[[',1))
    dat$celltype<-factor(dat$celltype, levels=c("Neurons","Astrocytes","Microglia","Endothelial","Oligodendrocytes",                                        "OPC","Fetal-quiescent","Fetal-replicating",
                                "Ex1","Ex2","Ex3","Ex4","Ex5","Ex6","Ex7","Ex8",
                                "In1","In2","In3","In4","In5","In6","In7","In8"))

    detach("package:dplyr", unload=TRUE)
    sedat<-summarySE(dat, measurevar="Expr", groupvars=c("celltype", "Group"))
    scale_this <- function(x) as.vector(scale(x))
    sedat$subtype<-ifelse(sedat$celltype %in% c(paste0("Ex",1:8), paste0("In",1:8)),"subtype","main")

    library(dplyr)
    scaledat = sedat %>% 
        group_by(Group, subtype) %>% 
        mutate(value=scale_this(Expr))

    return(scaledat)
}


scaledat<-target(neuron,glia,paste0(opt$Neuron,"_",opt$disease),paste0(opt$Glia,"_",opt$disease))

scaledat_main<-scaledat[scaledat$subtype=='main',]
scalemat1 = matrix(scaledat_main$value, nrow=2)
colnames(scalemat1) = unique(scaledat_main$celltype)
rownames(scalemat1) = scaledat_main$Group[1:2]

scaledat_sub<-scaledat[scaledat$subtype=='subtype',]
scalemat2 = matrix(scaledat_sub$value, nrow=2) 
colnames(scalemat2) = unique(scaledat_sub$celltype)
rownames(scalemat2) = scaledat_sub$Group[1:2]

col2plot = colorRampPalette(c("#053061","#2166AC","#4393C3","#92C5DE","#D1E5F0","#FFFFFF","#FDDBC7","#F4A582","#D6604D","#B2182B","#67001F"))

pdf(file=paste0(opt$outdir,"singlecell_exp_ind_",opt$disease,"_adjusted_005.pdf")) 

corrplot(scalemat1, is.corr=FALSE, method="color", cl.lim=c(-2.5,2.5), tl.col="black", col=col2plot(100)) 
corrplot(scalemat2, is.corr=FALSE, method="color", cl.lim=c(-2.5,2.5), tl.col="black", col=col2plot(100))
dev.off()
