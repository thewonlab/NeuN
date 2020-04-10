if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("Gviz")

library(Gviz)
library(S4Vectors)
library(stats4)
library(BiocGenerics)
library(parallel)
library(IRanges)
library(GenomicRanges)
library(GenomeInfoDb)
library(grid)
library(biomaRt)

FIRE<-read.table("FIRE_ANALYSIS_40KB.txt",header=T,sep=" ")
bm <- useMart(host="grch37.ensembl.org", biomart="ENSEMBL_MART_ENSEMBL", dataset="hsapiens_gene_ensembl")
NeuNpos<-FIRE[FIRE$NeuNpos_indicator==1,]
NeuNneg<-FIRE[FIRE$NeuNneg_indicator==1,]
pos_FIRE<-GRanges(seqnames=NeuNpos$chr, IRanges(start=NeuNpos$start, end=NeuNpos$end))
neg_FIRE<-GRanges(seqnames=NeuNneg$chr, IRanges(start=NeuNneg$start, end=NeuNneg$end))

itrack <- IdeogramTrack(genome="hg19", chromosome ="chr12")
gtrack <- GenomeAxisTrack(genome="hg19", chromosome ="chr12",
add53=TRUE,add35=TRUE,labelPos="below")

FIRE_pvalue<-FIRE_score<-GRanges(seqnames=FIRE$chr, IRanges(start=FIRE$start-20000, end=FIRE$end-20000))
values(FIRE_score) <- DataFrame(pos_score =FIRE[,4],neg_score =FIRE[,5])
values(FIRE_pvalue) <- DataFrame(pos_pvalue =FIRE[,6],neg_pvalue =FIRE[,7])
pos_atrack <- AnnotationTrack(lwd=1,stacking="dense",pos_FIRE, genome="hg19", chromosome="chr12", name="NeuNpos", rotate.title=360)
neg_atrack <- AnnotationTrack(lwd=1,stacking="dense",neg_FIRE, genome="hg19", chromosome="chr12", name="NeuNneg", rotate.title=360)
biomTrack <- BiomartGeneRegionTrack(genome="hg19",symbol="GRIN2B",
name="ENSEMBL", biomart=bm)
dTrack1 <- DataTrack(FIRE_score,genome="hg19", chromosome="chr12", type=c("s"), legend=TRUE,groups=c("NeuN+_score", "NeuN-_score"),col=c("bisque","blueviolet"),name="FIRE score")
dTrack2 <- DataTrack(FIRE_pvalue,genome="hg19", chromosome="chr12", type=c("s"), legend=TRUE,groups=c("NeuN+_-log10(pvalue)", "NeuN-_-log10(pvalue)"),col=c("brown1","darkcyan"),name="-log10(pvalue)")

pdf("FIRE_GRIN2B_total.pdf")
plotTracks(list(itrack,gtrack,biomTrack,pos_atrack,neg_atrack,dTrack1,dTrack2),from=13537337,to=14282012,
transcriptAnnotation="symbol",collapseTranscripts="meta",stackHeight=0.5,
background.title = "white",col.axis="black",col.title="black",cex.title=0.4,
cex.axis=0.3,cex.main=0.8)
dev.off()