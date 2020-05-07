### [1] Processing the data
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
library(Repitools)
setwd("/nas/longleaf/home/hubenxia/project/Loops/neuron_vs_glia")##  for loops dataset
load("geneAnno_allgenes.rda") # the saved file is geneAnno1
promoter<-read.table("gencode19_promoter.bed",sep="\t")
NeuNneg_TSS_interacting_region<-read.table("NeuNneg_TSS_interacting_region_fixed.bed", header=T,sep="\t")
NeuNneg_TSS_interacting_region_number<-nrow(NeuNneg_TSS_interacting_region)
NeuNneg_TSS_interacting_region_DF<-NeuNneg_TSS_interacting_region[,c(1,4,5)]
NeuNneg_loops_number<-nrow(unique(NeuNneg_TSS_interacting_region_DF))

NeuNpos_TSS_interacting_region = read.table("NeuNpos_TSS_interacting_region_fixed.bed", header=T,sep="\t")
NeuNpos_TSS_interacting_region_number<-nrow(NeuNpos_TSS_interacting_region)
NeuNpos_TSS_interacting_region_DF<-NeuNpos_TSS_interacting_region[,c(1,4,5)]
NeuNpos_loops_number<-nrow(unique(NeuNpos_TSS_interacting_region_DF))

TSS_ranges.NeuNpos<-GRanges(seqnames=NeuNpos_TSS_interacting_region[,1], IRanges(NeuNpos_TSS_interacting_region[,2], NeuNpos_TSS_interacting_region[,3]))
TSS_ranges.NeuNneg<-GRanges(seqnames=NeuNneg_TSS_interacting_region[,1], IRanges(NeuNneg_TSS_interacting_region[,2],NeuNneg_TSS_interacting_region[,3]))

mcols(TSS_ranges.NeuNpos)<-NeuNpos_TSS_interacting_region[ ,c(4:5)]
mcols(TSS_ranges.NeuNneg)<-NeuNneg_TSS_interacting_region[ ,c(4:5)]

promoteranges<-GRanges(seqnames=promoter[,1], IRanges(promoter[,2],promoter[,3]), gene=promoter[,6])

### [2] Overlap  for interactions with promoters
olap<-findOverlaps(TSS_ranges.NeuNpos, promoteranges)
TSS_pro.NeuNpos<-TSS_ranges.NeuNpos[queryHits(olap)]
mcols(TSS_pro.NeuNpos)<-cbind(mcols(TSS_pro.NeuNpos), mcols(promoteranges[subjectHits(olap)]))
TSS_pro.NeuNpos_DF<-annoGR2DF(TSS_pro.NeuNpos)
TSS_pro.NeuNpos_DF_fire<-TSS_pro.NeuNpos_DF[,c(1:3)]
TSS_pro.NeuNpos_DF_fire_number<-nrow(unique(TSS_pro.NeuNpos_DF_fire))


olap<-findOverlaps(TSS_ranges.NeuNneg, promoteranges)
TSS_pro.NeuNneg<-TSS_ranges.NeuNneg[queryHits(olap)]
mcols(TSS_pro.NeuNneg)<-cbind(mcols(TSS_pro.NeuNneg), mcols(promoteranges[subjectHits(olap)]))
TSS_pro.NeuNneg_DF<-annoGR2DF(TSS_pro.NeuNneg)
TSS_pro.NeuNneg_DF_fire<-TSS_pro.NeuNneg_DF[,c(1:3)]
TSS_pro.NeuNneg_DF_fire_number<-nrow(unique(TSS_pro.NeuNneg_DF_fire))


neg_TSS_DF <- unique(merge(TSS_pro.NeuNneg_DF,geneAnno1, by.x='gene',by.y='ensembl_gene_id')[,2:8])
pos_TSS_DF <- unique(merge(TSS_pro.NeuNpos_DF,geneAnno1, by.x='gene',by.y='ensembl_gene_id')[,2:8])
pos_TSS_DF <- pos_TSS_DF[pos_TSS_DF$hgnc_symbol != '',]
neg_TSS_DF <- neg_TSS_DF[neg_TSS_DF$hgnc_symbol != '',]
write.table(file="TSS_pro.NeuNneg_DF.txt",neg_TSS_DF,quote=F,row.names=F,sep="\t")
write.table(file="TSS_pro.NeuNpos_DF.txt",pos_TSS_DF,quote=F,row.names=F,sep="\t")

NeuNpos_TSS_interacting_nonpromoter_number<-NeuNpos_TSS_interacting_region_number-TSS_pro.NeuNpos_DF_fire_number
NeuNneg_TSS_interacting_nonpromoter_number<-NeuNneg_TSS_interacting_region_number-TSS_pro.NeuNneg_DF_fire_number

regions<-rep(c("TSS","nonTSS","TSS","nonTSS"),c(TSS_pro.NeuNpos_DF_fire_number,NeuNpos_TSS_interacting_nonpromoter_number,TSS_pro.NeuNneg_DF_fire_number,NeuNneg_TSS_interacting_nonpromoter_number))
sample_names<-rep(c("NeuNpos","NeuNpos","NeuNneg","NeuNneg"),c(TSS_pro.NeuNpos_DF_fire_number,NeuNpos_TSS_interacting_nonpromoter_number,TSS_pro.NeuNneg_DF_fire_number,NeuNneg_TSS_interacting_nonpromoter_number))
tss_loop<-data.frame(count=regions,type=sample_names)
#######################plot
pdf("geom_bar_tss_promoter.pdf", width=5, height=5)
ggplot(tss_loop,aes(x=type,fill=count))+
geom_bar(position ="fill")+
labs(title="proportion of TSS overlapping with Promoters \n all_tss_pos:187674 \n all_tss_neg:167551 \n tss_pos:15995 \n tss_neg:15898",y="proportion", x = "Hi-C samples")+
theme_bw() + theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank())
dev.off()

### [3] for gene symbol
TSS_NeuNposensg = unique(TSS_pro.NeuNpos$gene)
TSS_NeuNnegensg = unique(TSS_pro.NeuNneg$gene)
TSS_NeuNposgene = geneAnno1[geneAnno1$ensembl_gene_id %in% TSS_NeuNposensg, "hgnc_symbol"]
TSS_NeuNneggene = geneAnno1[geneAnno1$ensembl_gene_id %in% TSS_NeuNnegensg, "hgnc_symbol"]

TSS_NeuNposgene = TSS_NeuNposgene[TSS_NeuNposgene!=""]
TSS_NeuNneggene = TSS_NeuNneggene[TSS_NeuNneggene!=""]

write.table(file="TSS_NeuNposgene.txt",as.data.frame(TSS_NeuNposgene),quote=F,row.names=F)
write.table(file="TSS_NeuNneggene.txt",as.data.frame(TSS_NeuNneggene),quote=F,row.names=F)


###############for overlapping interactions with differential H3K27ac peaks
inter_ranges.NeuNpos = GRanges(seqnames=TSS_pro.NeuNpos_DF$chr, IRanges(TSS_pro.NeuNpos_DF$start_interaction, TSS_pro.NeuNpos_DF$end_interaction))
inter_ranges.NeuNneg = GRanges(seqnames=TSS_pro.NeuNneg_DF$chr, IRanges(TSS_pro.NeuNneg_DF$start_interaction, TSS_pro.NeuNneg_DF$end_interaction))

mcols(inter_ranges.NeuNpos)<-TSS_pro.NeuNpos_DF[ ,c("start","end","gene")]
mcols(inter_ranges.NeuNneg)<-TSS_pro.NeuNneg_DF[ ,c("start","end","gene")]

neuron<-read.table("neuron_H3K27ac_differential_peak.txt",sep="\t",header=T,quote="")
glia<-read.table("non-neuron_H3K27ac_differential_peak.txt",sep="\t",header=T,quote="")
############all differential peaks

neuronranges = GRanges(seqnames=neuron$Chrom, IRanges(neuron$Start, neuron$End))
gliaranges = GRanges(seqnames=glia$Chrom, IRanges(glia$Start, glia$End))

### [2] Overlap  for interactions with differential H3K27ac peaks
olap = findOverlaps(inter_ranges.NeuNpos, neuronranges)
inter_peak.NeuNpos = inter_ranges.NeuNpos[queryHits(olap)]
mcols(inter_peak.NeuNpos) = cbind(mcols(inter_peak.NeuNpos), annoGR2DF(neuronranges[subjectHits(olap)]))
inter_peak.NeuNpos_DF<-annoGR2DF(inter_peak.NeuNpos)


write.table(file="inter_peak.NeuNpos_DF_unique.txt",unique(inter_peak.NeuNpos_DF[,1:7]),sep="\t",quote=F,row.names=F)
write.table(file="inter_peak.NeuNpos_DF.txt",unique(inter_peak.NeuNpos_DF),sep="\t",quote=F,row.names=F)

######assigned differential peaks
inter_peak.NeuNpos_DF_peak<-inter_peak.NeuNpos_DF[,c(8:10)]
inter_peak.NeuNpos_DF_peak_number<-nrow(unique(inter_peak.NeuNpos_DF_peak))

neuron_enhancer_promoter<-nrow(unique(read.table("inter_peak.NeuNpos_DF_unique.txt",sep="\t",header=T)[,1:6]))
glia_enhancer_promoter<-nrow(unique(read.table("inter_peak.NeuNneg_DF_unique.txt",sep="\t",header=T)[,1:6]))

olap = findOverlaps(inter_ranges.NeuNneg, gliaranges)
inter_peak.NeuNneg = inter_ranges.NeuNneg[queryHits(olap)]
mcols(inter_peak.NeuNneg) = cbind(mcols(inter_peak.NeuNneg), annoGR2DF(gliaranges[subjectHits(olap)]))
inter_peak.NeuNneg_DF<-annoGR2DF(inter_peak.NeuNneg)

write.table(file="inter_peak.NeuNneg_DF_unique.txt",unique(inter_peak.NeuNneg_DF[,1:7]),sep="\t",quote=F,row.names=F)
write.table(file="inter_peak.NeuNneg_DF.txt",unique(inter_peak.NeuNneg_DF),sep="\t",quote=F,row.names=F)

pos<-read.table("inter_peak.NeuNpos_DF_unique.txt",sep="\t",header=T)
neg<-read.table("inter_peak.NeuNneg_DF_unique.txt",sep="\t",header=T)

pos<-merge(pos,geneAnno1,by.x='gene',by.y='ensembl_gene_id')[,c(2:4,6:8)]
neg<-merge(neg,geneAnno1,by.x='gene',by.y='ensembl_gene_id')[,c(2:4,6:8)]
pos<-pos[pos$hgnc_symbol!='',]
neg<-neg[neg$hgnc_symbol!='',]
write.table(file="TSS_pro.NeuNpos_DF_name.bed",pos,sep="\t",quote=F,row.names=F,col.names=F)
write.table(file="TSS_pro.NeuNneg_DF_name.bed",neg,sep="\t",quote=F,row.names=F,col.names=F)

inter_peak.NeuNneg_DF_peak<-inter_peak.NeuNneg_DF[,c(8:10)]
inter_peak.NeuNneg_DF_peak_number<-nrow(unique(inter_peak.NeuNneg_DF_peak))

inter_peak.NeuNpos_DF_fire<-inter_peak.NeuNpos_DF[,c(1:3)]
inter_peak.NeuNpos_DF_fire_number<-nrow(unique(inter_peak.NeuNpos_DF_fire))
inter_peak.NeuNneg_DF_fire<-inter_peak.NeuNneg_DF[,c(1:3)]
inter_peak.NeuNneg_DF_fire_number<-nrow(unique(inter_peak.NeuNneg_DF_fire))

NeuNpos_interacting_nonpeak_number<-NeuNpos_TSS_interacting_region_number-inter_peak.NeuNpos_DF_fire_number
NeuNneg_interacting_nonpeak_number<-NeuNneg_TSS_interacting_region_number-inter_peak.NeuNneg_DF_fire_number

regions<-rep(c("Loops","nonLoops","Loops","nonLoops"),c(inter_peak.NeuNpos_DF_fire_number,NeuNpos_interacting_nonpeak_number,inter_peak.NeuNneg_DF_fire_number,NeuNneg_interacting_nonpeak_number))
sample_names<-rep(c("NeuNpos","NeuNpos","NeuNneg","NeuNneg"),c(inter_peak.NeuNpos_DF_fire_number,NeuNpos_interacting_nonpeak_number,inter_peak.NeuNneg_DF_fire_number,NeuNneg_interacting_nonpeak_number))
loop_peak<-data.frame(count=regions,type=sample_names)
#######################plot
pdf("geom_bar_loops_H3K27ac.pdf", width=5, height=5)
ggplot(loop_peak,aes(x=type,fill=count))+
geom_bar(position ="fill")+
labs(title="proportion of Loops overlapping with peaks \n all_Loops_pos:187674 \n all_Loops_neg:167551 \n Loops_pos:16370 \n Loops_neg:19418",y="proportion", x = "Hi-C samples")+
theme_bw() + theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank())
dev.off()

### [3] for gene symbol
inter_peak_NeuNposensg = unique(inter_peak.NeuNpos$gene)
inter_peak_NeuNnegensg = unique(inter_peak.NeuNneg$gene)
inter_peak_NeuNposgene = geneAnno1[geneAnno1$ensembl_gene_id %in% inter_peak_NeuNposensg, "hgnc_symbol"]
inter_peak_NeuNneggene = geneAnno1[geneAnno1$ensembl_gene_id %in% inter_peak_NeuNnegensg, "hgnc_symbol"]

inter_peak_NeuNposgene = inter_peak_NeuNposgene[inter_peak_NeuNposgene!=""]
inter_peak_NeuNneggene = inter_peak_NeuNneggene[inter_peak_NeuNneggene!=""]

write.table(file="inter_peak_NeuNposgene.txt",as.data.frame(unique(inter_peak_NeuNposgene)),quote=F,row.names=F)
write.table(file="inter_peak_NeuNneggene.txt",as.data.frame(unique(inter_peak_NeuNneggene)),quote=F,row.names=F)

###################for heatmap
targetgene = list(unique(inter_peak_NeuNposgene),unique(inter_peak_NeuNneggene))
diseasename = c("NeuNpos", "NeuNneg")
cellexp<-read.table("Capstone4_singlecellexpression_adultbrain.txt", header=T)
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

dat$cell = gsub("Fetal.quiescent","Fetal-quiescent",dat$cell)
dat$cell = gsub("Fetal.replicating","Fetal-replicating",dat$cell)

dat$celltype = unlist(lapply(strsplit(dat$cell, split="[.]"),'[[',1))

dat = dat[!(dat$celltype %in% c("OPC","Fetal-quiescent","Fetal-replicating",
                              "Ex1","Ex2","Ex3","Ex4","Ex5","Ex6","Ex7","Ex8",
                              "In1","In2","In3","In4","In5","In6","In7","In8")), ]

dat$celltype = factor(dat$celltype, levels=c("Neurons","Astrocytes","Microglia","Endothelial","Oligodendrocytes"))

pdf("Loops_matched_gene_cell-type_expression.pdf", width=5, height=5)
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
pdf(file="singlecell_exp_NeuNposVSneg_Loops.pdf", width=10, height=5)
corrplot(scalemat, is.corr=FALSE, method="color", tl.col="black", col=col2plot(100), cl.lim=c(-2.0,2.0)) # 
dev.off()
