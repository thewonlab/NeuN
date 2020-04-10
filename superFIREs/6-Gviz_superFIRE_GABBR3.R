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
NeuNpos<-read.table("super_FIRE_call_NeuNpos.txt",sep="\t",header=T)
NeuNneg<-read.table("super_FIRE_call_NeuNneg.txt",sep="\t",header=T)
pos_FIRE<-GRanges(seqnames=NeuNpos$chr, IRanges(start=NeuNpos$start, end=NeuNpos$end))
neg_FIRE<-GRanges(seqnames=NeuNneg$chr, IRanges(start=NeuNneg$start, end=NeuNneg$end))

itrack <- IdeogramTrack(genome="hg19", chromosome ="chr15")
gtrack <- GenomeAxisTrack(genome="hg19", chromosome ="chr15",
add53=TRUE,add35=TRUE,labelPos="below")

FIRE_pvalue<-GRanges(seqnames=FIRE$chr, IRanges(start=FIRE$start-20000, end=FIRE$end-20000))
values(FIRE_pvalue) <- DataFrame(pos_pvalue =FIRE[,6],neg_pvalue =FIRE[,7])
pos_atrack <- AnnotationTrack(lwd=1,stacking="dense",pos_FIRE, genome="hg19", chromosome="chr15", name="NeuNpos", rotate.title=360)
neg_atrack <- AnnotationTrack(lwd=1,stacking="dense",neg_FIRE, genome="hg19", chromosome="chr15", name="NeuNneg", rotate.title=360)
biomTrack <- BiomartGeneRegionTrack(genome="hg19",symbol="GABRB3",
name="ENSEMBL", biomart=bm)
dTrack <- DataTrack(FIRE_pvalue,genome="hg19", chromosome="chr15", type=c("s"), legend=TRUE,groups=c("NeuN+_-log10(pvalue)", "NeuN-_-log10(pvalue)"),col=c("brown1","darkcyan"),name="-log10(pvalue)")

pdf("superFIRE_GABRB3.pdf")
plotTracks(list(itrack,gtrack,biomTrack,pos_atrack,neg_atrack,dTrack),from=26500000,to=28000000,
transcriptAnnotation="symbol",collapseTranscripts="meta",stackHeight=0.5,
background.title = "white",col.axis="black",col.title="black",cex.title=0.4,
cex.axis=0.3,cex.main=0.8)
dev.off()