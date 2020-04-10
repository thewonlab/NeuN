options(stringsAsFactors=F)
library(stats4)
library(BiocGenerics)
library(parallel)
library(S4Vectors)
library(IRanges)
library(GenomeInfoDb)
library(GenomicRanges)
library(ggplot2)
library(reshape2)
library(dplyr)

setwd("/nas/longleaf/home/hubenxia/project/Loops/GLU_GABA/")
GLUdifferential<-read.table("GLUdiff_promoter_enhancer_loops_peak_gene.txt",header=F)[,1]
GABAdifferential<-read.table("GABAdiff_promoter_enhancer_loops_peak_gene.txt",header=F)[,1]
GLUspecificity<-read.table("GLU_promoter_enhancer_loops_cell_type_peak_gene.txt",header=F)[,1]
GABAspecificity<-read.table("GABA_promoter_enhancer_loops_cell_type_peak_gene.txt",header=F)[,1]
GLULoops<-unique(c(GLUdifferential,GLUspecificity))
GABALoops<-unique(c(GABAdifferential,GABAspecificity))
GLUpromoter<-read.table("GLUdiff_promoter_gene.txt",header=F)[,1]
GABApromoter<-read.table("GABAdiff_promoter_gene.txt",header=F)[,1]
targetgene<-list(GLUdifferential,GABAdifferential,GLUspecificity,GABAspecificity,GLULoops,GABALoops,GLUpromoter,GABApromoter)
diseasename<-c('GLUdifferential','GABAdifferential','GLUspecificity','GABAspecificity','GLULoops','GABALoops','GLUpromoter','GABApromoter')

cellexp<-read.table("/nas/longleaf/home/hubenxia/project/Loops/neuron_vs_glia/Capstone4_singlecellexpression_adultbrain.txt", header=T)
rownames(cellexp) = cellexp[,1]
cellexp = cellexp[,-1]
datExpr = scale(cellexp,center=T, scale=F) 
exprdat = vector(mode="list", length=length(targetgene))

for(i in 1:length(targetgene)){
  exprdat[[i]] = apply(datExpr[match(targetgene[[i]], rownames(datExpr)),],2,mean,na.rm=T)
}

dat = c()

for(i in 1:length(targetgene)){
  datframe = data.frame(Group=diseasename[i], cell=names(exprdat[[i]]), Expr=exprdat[[i]])
  dat = rbind(dat, datframe)
}

dat$cell = gsub("fetal.quiescent","fetal-quiescent",dat$cell)
dat$cell = gsub("fetal.replicating","fetal-replicating",dat$cell)

dat$celltype = unlist(lapply(strsplit(dat$cell, split="[.]"),'[[',1))

dat = (dat[(dat$celltype %in% c("Ex1","Ex2","Ex3","Ex4","Ex5","Ex6","Ex7","Ex8","In1","In2","In3","In4","In5","In6","In7","In8")),])

dat$celltype = factor(dat$celltype, levels=c("Ex1","Ex2","Ex3","Ex4","Ex5","Ex6","Ex7","Ex8","In1","In2","In3","In4","In5","In6","In7","In8"))

pdf("GLU_GABA_gene_cell-type_expression.pdf", width=5, height=5)
ggplot(dat,aes(x=celltype, y=Expr, fill=Group)) + geom_hline(yintercept=0, colour="red") + #scale_fill_manual(values=c("#99CA3C","#75B3E1")) + 
  ylab("Normalized expression") +  geom_boxplot(outlier.size = NA) + theme(axis.text.x=element_text(angle = 90, hjust=1)) +
  ggtitle("Expression in each cell type: Loops") 
dev.off()

library(corrplot)
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  library(plyr)

  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) {
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
  
  # Rename the "mean" column    
  datac <- rename(datac, c("mean" = measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}
detach("package:dplyr", unload=TRUE)

sedat = summarySE(dat, measurevar="Expr", groupvars=c("celltype", "Group"))

scale_this <- function(x) as.vector(scale(x))
library(dplyr)
scaledat = sedat %>% 
  group_by(Group) %>% 
  mutate(value=scale_this(Expr)) 


scalemat = matrix(scaledat$value, nrow=length(diseasename)) 
colnames(scalemat) = unique(scaledat$celltype)
rownames(scalemat) = scaledat$Group[1:length(diseasename)]

col2plot = colorRampPalette(c("#053061","#2166AC","#4393C3","#92C5DE","#D1E5F0","#FFFFFF","#FDDBC7","#F4A582","#D6604D","#B2182B","#67001F"))
pdf(file="singlecell_exp_GLU_GABA.pdf", width=10, height=5)
corrplot(scalemat, is.corr=FALSE, method="color", tl.col="black", col=col2plot(100), cl.lim=c(-2,2)) # 
dev.off()

col2plot = colorRampPalette(c("#053061","#2166AC","#4393C3","#92C5DE","#D1E5F0","#FFFFFF","#FDDBC7","#F4A582","#D6604D","#B2182B","#67001F"))
pdf(file="singlecell_exp_GLU_GABA_EX.pdf", width=10, height=5)
corrplot(scalemat[,1:8], is.corr=FALSE, method="color", tl.col="black", col=col2plot(100), cl.lim=c(-2,2)) # 
dev.off()

col2plot = colorRampPalette(c("#053061","#2166AC","#4393C3","#92C5DE","#D1E5F0","#FFFFFF","#FDDBC7","#F4A582","#D6604D","#B2182B","#67001F"))
pdf(file="singlecell_exp_GLU_GABA_IN.pdf", width=10, height=5)
corrplot(scalemat[,9:16], is.corr=FALSE, method="color", tl.col="black", col=col2plot(100), cl.lim=c(-2,2)) # 
dev.off()