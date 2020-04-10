### [1] Get the distribution of the distance between interaction
options(stringsAsFactors=F)
library(GenomicRanges)
library(ggplot2)
setwd("/proj/hyejunglab/NeuN/Loop")
tssneg = read.table("NeuNneg_TSS_interacting_region_fixed.bed",header=T)
tsspos = read.table("NeuNpos_TSS_interacting_region_fixed.bed",header=T)
tssneg$distance = abs(tssneg[,2]-tssneg[,4])
tsspos$distance = abs(tsspos[,2]-tsspos[,4])

pdf("TSS_NeuN_distance_distribution.pdf", width=10, height=5)
par(mfrow=c(1,2))
hist(tssneg$distance, xlim=c(0,1000000), main="NeuN- TSS",xlab="Distance (bp)",col="midnightblue",border="white",ylab="Frequency")
abline(h=0,v=median(tssneg$distance),col="red")  #350000
sum(tssneg$distance<median(tssneg$distance))/length(tssneg$distance) # 0.4903403
hist(tsspos$distance, xlim=c(0,1000000), main="NeuN+ TSS",xlab="Distance (bp)",col="midnightblue",border="white",ylab="Frequency")
abline(h=0,v=median(tsspos$distance),col="red")   #320000
sum(tsspos$distance<median(tsspos$distance))/length(tsspos$distance) # 0.4993446
dev.off()

### [2] Functional enrichment of TSS interacting regions

chrompeak = "/proj/hyejunglab/epigenetics/Enhancer/chromHMM/core15/E073_15_coreMarks_dense.bed"
genelength = 2897310462 # http://genomewiki.ucsc.edu/index.php/Hg19_Genome_size_statistics
hicnegranges = GRanges(seqnames=tssneg[,1],ranges=IRanges(as.numeric(tssneg[,2]),as.numeric(tssneg[,3])))
hicposranges = GRanges(seqnames=tsspos[,1],ranges=IRanges(as.numeric(tsspos[,2]),as.numeric(tsspos[,3])))

##### ChromHMM 
chrom = read.table(chrompeak, skip=1); 
chrommat = data.frame("chr"=chrom[,1], "start"=chrom[,2], "end"=chrom[,3], "mark"=chrom[,4])
chrommat1 = chrommat[!(chrommat$chr %in% c("chrM", "chrY")), ]
chrommat = chrommat1
chromranges = GRanges(seqnames=chrommat$chr,ranges=IRanges(as.numeric(chrommat$start),as.numeric(chrommat$end)), 
                      mark = chrommat$mark);
chrommark<-sort(unique(chrommat$mark))
chrommark<-chrommark[c(7:15, 1:6)]

overlapping<-function(hicranges,chromranges,chrommark,filename){
    olap<-findOverlaps(hicranges,chromranges);
    hicvalues<-hicranges[queryHits(olap)];
    mcols(hicvalues)<-cbind(mcols(hicranges[queryHits(olap)]), mcols(chromranges[subjectHits(olap)]))
    pval<-c()
    enrichment<-c()
    for(i in 1:length(chrommark)){
        juicemark<-unique(hicvalues[hicvalues$mark==chrommark[i]])
        hisrangemark<-chromranges[chromranges$mark==chrommark[i]]
        p.binom<-sum(as.numeric(ranges(hisrangemark)@width))/genelength ## bp # of histone marks/total gene length  
        n.binom<-length(hicranges) # CAVIAR SNP #
        s.binom<-length(juicemark) # overlapping # of CAVIAR SNPs in DHS 
        pval<-c(pval, pbinom(s.binom-1,n.binom,p.binom,lower.tail=FALSE,log.p=TRUE) * -log10(exp(1)))
        numerator<-sum(width(reduce(juicemark)))/genelength;
        denominator<-(sum(as.numeric(width(reduce(hisrangemark))))/genelength) * (sum(width(reduce(hicranges)))/genelength);
        foldenrichment<-numerator/denominator;
        foldenrichment[which(foldenrichment==0)]<-NA;
        enrichment<-c(enrichment, foldenrichment)
    }
    names(pval)<-names(enrichment)<-chrommark
    pvalchrom<-pval
    pvalhisg<-data.frame("EpigeneticMarks"=names(pval), "Pval"=pval)
    pdf(file=paste0("histonemark_Pval_",filename,".pdf"), width=5, height=4)
    marker<-c('TssA','TssAFlnk','TxFlnk','Tx','TxWk','EnhG','Enh','ZNF/Rpts','Het',
            'TssBiv','BivFlnk','EnhBiv','ReprPC','ReprPCWk','Quies')
    label<-paste0(c(1:15),'_',marker)
    print(label)
    print(head(pvalhisg))
    pvalhisg$EpigeneticMarks<-factor(pvalhisg$EpigeneticMarks,level=label)
    barpl<-ggplot(pvalhisg, aes(x=EpigeneticMarks,y=Pval)) + 
    geom_bar(stat="identity",position="dodge",fill="midnightblue") + theme_bw() + 
    theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(), axis.text.x=element_text(angle=90,hjust=1)) + ylab("-log10 P-value") + xlab("")
    print(barpl)
    dev.off()

    pvalchromg<-data.frame("EpigeneticMarks"=names(pvalchrom), "Pval"=pvalchrom)
    positions <- names(pvalchrom)
    pvalchromg$EpigeneticMarks<-factor(pvalchromg$EpigeneticMarks,level=label)
    pdf(file=paste0("epimark_Pval_",filename,".pdf"), width=6, height=4)
    barp2<-ggplot(pvalchromg, aes(x=factor(EpigeneticMarks,level=label),y=Pval)) + 
    geom_bar(stat="identity",position="dodge",fill="midnightblue") + theme_bw() + scale_x_discrete(limits = positions) + 
    theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(), axis.text.x=element_text(angle=90,hjust=1)) + ylab("-log10 P-value") + xlab("")
    print(barp2)
    dev.off()

    pvalchrome<-data.frame("EpigeneticMarks"=names(enrichment), "Fold_enrichment"=enrichment)
    positions<-names(pvalchrom)
    pvalchrome$EpigeneticMarks<-factor(pvalchrome$EpigeneticMarks,level=label)
    pdf(file=paste0("epimark_enrichment_",filename,".pdf"), width=6, height=4)
    barp3<-ggplot(pvalchrome, aes(x=factor(EpigeneticMarks,level=label),y=Fold_enrichment))+ 
    geom_bar(stat="identity",position="dodge",fill="midnightblue") + theme_bw() + scale_x_discrete(limits = positions) + 
    theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(), axis.text.x=element_text(angle=90,hjust=1)) + ylab("Fold enrichment") + xlab("")
    print(barp3)
    dev.off()
}


overlapping(hicnegranges,chromranges,chrommark,'NeuNneg')
overlapping(hicposranges,chromranges,chrommark,'NeuNpos')















olap = findOverlaps(hicnegranges,chromranges);
hicvalues = hicnegranges[queryHits(olap)];
mcols(hicvalues) = cbind(mcols(hicnegranges[queryHits(olap)]), mcols(chromranges[subjectHits(olap)]))

chrommark = sort(unique(chrommat$mark))
chrommark = chrommark[c(7:15, 1:6)]
# chrommark = chrommark[c(11,18:25,1:10,12:17)] # for 25 marks

pval = c()
enrichment = c()
for(i in 1:length(chrommark)){
  juicemark = unique(hicvalues[hicvalues$mark==chrommark[i]])
  hisrangemark = chromranges[chromranges$mark==chrommark[i]]
  p.binom = sum(as.numeric(ranges(hisrangemark)@width))/genelength ## bp # of histone marks/total gene length  
  n.binom = length(hicnegranges) # CAVIAR SNP #
  s.binom = length(juicemark) # overlapping # of CAVIAR SNPs in DHS 
  
  pval = c(pval, pbinom(s.binom-1,n.binom,p.binom,lower.tail=FALSE,log.p=TRUE) * -log10(exp(1)))
  
  numerator = sum(width(reduce(juicemark)))/genelength;
  denominator = (sum(as.numeric(width(reduce(hisrangemark))))/genelength) * (sum(width(reduce(hicnegranges)))/genelength);
  foldenrichment = numerator/denominator;
  foldenrichment[which(foldenrichment==0)] = NA;
  enrichment = c(enrichment, foldenrichment)
}

names(pval) = names(enrichment) = chrommark

pvalchrom = pval

pvalhisg = data.frame("EpigeneticMarks"=names(pval), "Pval"=pval)
#positions = c("DHS", "H3K4me3", "H3K4me1", "H3K36me3", "H3K9me3", "H3K27me3")
pdf(file="histonemark_Pval_fetalbrain_NeuNbulk.pdf", width=5, height=4)
barpl = ggplot(pvalhisg, aes(x=EpigeneticMarks,y=Pval))
barpl + geom_bar(stat="identity",position="dodge",fill="midnightblue") + theme_bw() + 
  theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(), axis.text.x=element_text(angle=90,hjust=1)) + ylab("-log10 P-value") + xlab("")
dev.off()

pvalchromg = data.frame("EpigeneticMarks"=names(pvalchrom), "Pval"=pvalchrom)
positions <- names(pvalchrom)
pdf(file="epimark_Pval_fetalbrain_NeuNbulk.pdf", width=6, height=4)
barpl = ggplot(pvalchromg, aes(x=EpigeneticMarks,y=Pval))
barpl + geom_bar(stat="identity",position="dodge",fill="midnightblue") + theme_bw() + scale_x_discrete(limits = positions) + 
  theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(), axis.text.x=element_text(angle=90,hjust=1)) + ylab("-log10 P-value") + xlab("")
dev.off()

pvalchrome = data.frame("EpigeneticMarks"=names(enrichment), "Fold_enrichment"=enrichment)
positions <- names(pvalchrom)
pdf(file="epimark_enrichment_fetalbrain_NeuNbulk.pdf", width=6, height=4)
barpl = ggplot(pvalchrome, aes(x=EpigeneticMarks,y=Fold_enrichment))
barpl + geom_bar(stat="identity",position="dodge",fill="midnightblue") + theme_bw() + scale_x_discrete(limits = positions) + 
  theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(), axis.text.x=element_text(angle=90,hjust=1)) + ylab("Fold enrichment") + xlab("")
dev.off()