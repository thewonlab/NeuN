options(stringsAsFactors=F)
library(ggplot2)
library(GenomicRanges)
library(plyr)
setwd("/proj/hyejunglab/NeuN/Loop")
load("geneAnno_allgenes.rda")
promoter<-read.table("/proj/hyejunglab/chr/gencode19_promoter.bed")
promoter<-unique(promoter[,c(1,2,3,6)])
promoteranges = GRanges(seqnames=promoter[,1], IRanges(promoter[,2], promoter[,3]), gene=promoter[,4])

NeuNneg_TSS_interacting_region<-read.table("NeuNneg_TSS_interacting_region_fixed.bed", header=T,sep="\t")
NeuNpos_TSS_interacting_region = read.table("NeuNpos_TSS_interacting_region_fixed.bed", header=T,sep="\t")

peak<-read.table("H3K27ac.txt",sep="\t",header=T,quote="")
nondiffpeak<-peak[peak[,5]>0.05,]
neuronpeak<-read.table("neuron_H3K27ac_differential_peak.txt",header=T,sep="\t",quote="")
gliapeak<-read.table("non-neuron_H3K27ac_differential_peak.txt",header=T,sep="\t",quote="")
pospeak<-rbind(nondiffpeak,neuronpeak)[,c(2:4)]
negpeak<-rbind(nondiffpeak,gliapeak)[,c(2:4)]

condata<-function(chromAB){
    colnames(chromAB) =c("chr", "start", "end")
    chromAB = chromAB[!(chromAB$chr %in% c("chrY", "chrM")),]
    abranges = GRanges(seqnames=chromAB$chr, IRanges(chromAB$start, chromAB$end))
    return(abranges)
}

poschrom =condata(pospeak)
negchrom =condata(negpeak)


grange<-function(abint,promoteranges,abranges,filename){
    abregion = GRanges(seqnames=abint[,1], IRanges(abint[,2], abint[,2]+10000), chrint=abint[,4])
    #abdist = unique(data.frame(tss=paste(abint[,1],abint[,4],sep=":"), enh=paste(abint[,2],abint[,4],sep=":")))
    olap = findOverlaps(abregion, promoteranges)
    abpromoter = abregion[queryHits(olap)]
    mcols(abpromoter) = cbind(mcols(abregion[queryHits(olap)]),mcols(promoteranges[subjectHits(olap)]))
    abenhancer = GRanges(seqnames=seqnames(abpromoter), IRanges(abpromoter$chrint, abpromoter$chrint+10000), tss=start(abpromoter), gene=abpromoter$gene)
    olap = findOverlaps(abenhancer, abranges)
    abchromranges = abenhancer[queryHits(olap)]
    mcols(abchromranges) = cbind(mcols(abenhancer[queryHits(olap)]),mcols(abranges[subjectHits(olap)]))
    # To get the individual enhancer-promoter connections
    mcols2add = data.frame(mcols(abenhancer[queryHits(olap)]), capstone4.enh.start=start(abranges[subjectHits(olap)]), capstone4.enh.end=end(abranges[subjectHits(olap)])) 
    mcols(abchromranges) = mcols2add
    abEPdat = unique(data.frame(chr=seqnames(abchromranges), tss=abchromranges$tss, gene=abchromranges$gene, enh.start=abchromranges$capstone4.enh.start, enh.end=abchromranges$capstone4.enh.end))
    abEPdat$tss.start = abEPdat$tss
    abEPdat$tss.end = abEPdat$tss + 10000
    abEPdat = abEPdat[,-2]
    save(abEPdat, file=paste0("enhancer-promoter_interactions_",filename,".rda"))

    abactranges = abchromranges
    abactgenes = unique(abactranges$gene)
    abact = unique(data.frame(enh=paste(seqnames(abactranges),start(abactranges),sep=":"),tss=paste(seqnames(abactranges),abactranges$tss,sep=":")))
    abactcount = count(abact, "tss")
    abactcount[abactcount$freq>=10, "freq"] = 10

    dat = data.frame(table(abactcount$freq))

    enhnumplot = ggplot(dat, aes(x=Var1, y=Freq, group=1)) + theme_classic() + 
    geom_point(size=3, col="darkolivegreen3") + geom_line() + # geom_segment(aes(x=Var1, xend=Var1, y=0, yend=Freq)) + 
    ggtitle("") + xlab("# of enhancers interacting with promoters") + ylab("Frequency")
  
    pdf(paste0("enhancer-promoter_counts_",filename,".pdf"), height=5, width=6)
    p<-hist(abactcount$freq, breaks=0:10, xlab="# of enhancers interacting with promoters", ylab="Frequency", main="", col="midnightblue", border="white")
    print(enhnumplot)
    print(p)
    dev.off()
    return(abactranges)
}
posgene=grange(NeuNpos_TSS_interacting_region,promoteranges,poschrom,'NeuNpos')
neggene=grange(NeuNneg_TSS_interacting_region,promoteranges,negchrom,'NeuNneg')



## Will output this for future use: to share with other people
#abpromdat = data.frame(chr=seqnames(abpromoter), enhancer.start=as.integer(abpromoter$chrint), enhancer.end=as.integer(abpromoter$chrint+10000), promoter.start=as.integer(start(abpromoter)), gene=abpromoter$gene)
#abpromdat = unique(abpromdat)
#write.table(abpromdat, file="/proj/hyejunglab/TSS/TSS_interacting_region_AB_genemapped.bed", sep="\t", quote=F, row.names=F, col.names=T)


### [2] Clean the expression dataset from Capstone 4 
datExpr = read.table("/proj/hyejunglab/expression/capstone4/MERGED_Bipseq_BrainGVEX_Brainspan_CMC_CMC_HBCC_GTEx_LIBD_UCLA_Yale_QUANTILE_normalized_min_50.expression.txt", header=T)
rownames(datExpr) = datExpr[,1]
datExpr = datExpr[,-1]
gtexdatexpr = datExpr[,grep("GTEX",colnames(datExpr))]
gtexnames = colnames(gtexdatexpr)
gtexnames = paste0(unlist(lapply(strsplit(gtexnames, split="[.]"), '[[', 1)), "-", unlist(lapply(strsplit(gtexnames, split="[.]"), '[[', 2)))
colnames(datExpr)[grep("GTEX",colnames(datExpr))] <- gtexnames
indid = colnames(datExpr) # individual IDs
rownames(datExpr) = unlist(lapply(strsplit(rownames(datExpr), split="[.]"), '[[', 1))


### [3] For enhancer-promoter interactions, let's subgroup the genes based on the # of enhancers they interact with

expr<-function(datExpr,abactranges,filename){
    datExpr = scale(log2(datExpr+1),center = T, scale=FALSE) + 1
    abactbase = unique(data.frame(gene=abactranges$gene,enh=start(abactranges)))
    abactcount = count(abactbase, "gene")
    N=10
    abactcount[abactcount$freq>=N, "freq"] = N
    exprdat = vector(mode="list", length=N)
    for(i in 1:N){
        targetgene = abactcount[abactcount$freq==i, "gene"]
        #print(length(targetgene))
        exprdat[[i]] = apply(datExpr[match(targetgene, rownames(datExpr)),],2,mean,na.rm=T)
        #print(range(exprdat[[i]]))
    }
    dat = c()
    for(i in 1:N){
        datframe = data.frame(Group=i,Expr=exprdat[[i]])
        dat = rbind(dat, datframe)
    }
    print(head(dat))
    p.brainExp = 
    ggplot(dat,aes(x=as.factor(Group), y=Expr)) + 
    ylab("Normalized expression") +  geom_boxplot(outlier.size = NA, fill="midnightblue")  +
    ggtitle(paste0(filename," Expression")) + xlab("# of enhancers") 
    pdf(paste0("enhancer_expression_",filename,".pdf"), width=5, height=5)
    print(p.brainExp)
    dev.off()
    pdf(paste0(filename,"_loops_enhancer_number_expression.pdf"))
    p <- boxplot(Expr ~ Group,dat)
    new_data<-aggregate(Expr ~ Group,data=dat, FUN=mean)
    fit<-lm(Expr ~ Group,data=new_data)
    abline(fit,col='red')
    print(summary(fit))
    print(p)
    dev.off()
}
expr(datExpr,posgene,'NeuNpos')
#Call:
lm(formula = Expr ~ Group, data = new_data)

Residuals:
     Min       1Q   Median       3Q      Max 
-0.25123 -0.03101  0.01011  0.06765  0.19644 

Coefficients:
            Estimate Std. Error t value Pr(>|t|)    
(Intercept)  1.35184    0.09212  14.674 4.57e-07 ***
Group        0.08360    0.01485   5.631 0.000492 ***
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Residual standard error: 0.1349 on 8 degrees of freedom
Multiple R-squared:  0.7985,	Adjusted R-squared:  0.7734 
F-statistic: 31.71 on 1 and 8 DF,  p-value: 0.0004921



expr(datExpr,neggene,'NeuNneg')

#Call:
lm(formula = Expr ~ Group, data = new_data)

Residuals:
     Min       1Q   Median       3Q      Max 
-0.27200 -0.03443  0.01483  0.07132  0.12929 

Coefficients:
            Estimate Std. Error t value Pr(>|t|)    
(Intercept)  1.28343    0.08516   15.07 3.72e-07 ***
Group        0.06917    0.01372    5.04    0.001 ** 
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Residual standard error: 0.1247 on 8 degrees of freedom
Multiple R-squared:  0.7605,	Adjusted R-squared:  0.7305 
F-statistic:  25.4 on 1 and 8 DF,  p-value: 0.001002