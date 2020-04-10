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
library(rtracklayer)

setwd("/nas/longleaf/home/hubenxia/project/Loops/GLU_GABA/")
#bigwig directoy: /proj/hyejunglab/epigenetics/NeuNsort/EpiMap/bigwig
##Hi-C data
hic_GABA<-read.table("GABAdiff_promoter_enhancer_loops_peak.txt",header=F,sep="\t")
hic_GLU<-read.table("GLUdiff_promoter_enhancer_loops_peak.txt",header=F,sep="\t")

enhrange<-function(input_file){
    enhancer_range<-GRanges(seqnames=input_file[,1],IRanges(input_file[,4],input_file[,5]),tss=input_file[,2],gene=input_file[,6])
    return(enhancer_range)
}

prorange<-function(input_file){
    promoter_range<-GRanges(seqnames=input_file[,1],IRanges(input_file[,2],input_file[,3]),tss=input_file[,2],gene=input_file[,6])
    return(promoter_range)
}

hicenhancer_GABA<-enhrange(hic_GABA)
hicpromoter_GABA<-prorange(hic_GABA)
hicenhancer_GLU<-enhrange(hic_GLU)
hicpromoter_GLU<-prorange(hic_GLU)

##Generate other required files
tad_GABA<-read.table("/proj/hyejunglab/NeuN/NeuNpos/chrmat/HMMdomain.bed")
tad_GLU<-read.table("/proj/hyejunglab/NeuN/NeuNpos/chrmat/HMMdomain.bed")
tad_range<-function(tad){
    tadranges<-GRanges(seqnames=tad[,1], IRanges(tad[,2], tad[,3]))
    return(tadranges)
}
tad_GABA_ran<-tad_range(tad_GABA)
tad_GLU_ran<-tad_range(tad_GLU)

############all differential peaks
peak <- read.table("/proj/hyejunglab/NeuN/ChIPseq/diffbind/H3K27ac3_differential_peaks.txt",sep="\t",header=T)
GLU_peak <- peak[peak$Fold >= 1 & peak$FDR<=0.01,][,c(1:3)]
GABA_peak <- peak[peak$Fold <= -1 & peak$FDR<=0.01,][,c(1:3)]

GLUranges<-tad_range(GLU_peak)
GABAranges<-tad_range(GABA_peak)


load("/proj/hyejunglab/chr/geneAnno_allgenes.rda")
load("/proj/hyejunglab/chr/itrackchromosome.rda")
mart <- useMart(biomart="ENSEMBL_MART_ENSEMBL",dataset="hsapiens_gene_ensembl",host="feb2014.archive.ensembl.org")

#define Neuronal genes vs. glial genes
genelist<-c("GRIK4","GLUL","GRIA1","SHANK2","PRDM2","GABRB1") 

Track_new<-function(enh_range,pro_range,hic_name,chr_name){
    chrnum<-chr_name # chromosome number
    thisHiC.Enh<-enh_range
    thisHiC.TSS<-pro_range
    thisHiC.Enh<-thisHiC.Enh[thisHiC.Enh$gene==genelist[1]]  ##we need to ensure thatthisHiC.Enh is not NULL
    thisHiC.TSS<-thisHiC.TSS[thisHiC.TSS$gene==genelist[1]]  ##we need to ensure thatthisHiC.TSS is not NULL
    if(length(thisHiC.Enh) & length(thisHiC.TSS)){
        combined<-paste(start(thisHiC.Enh),end(thisHiC.Enh),start(thisHiC.TSS),end(thisHiC.TSS))
        dupind<-which(!duplicated(combined))
        thisHiC.Enh<-thisHiC.Enh[dupind]
        thisHiC.TSS<-thisHiC.TSS[dupind]

    ##Combine these in a way that sashimi will plot them as interactions
        thisHiC<-thisHiC.Enh
        mcols(thisHiC)<-NULL
        for (j in 1:length(thisHiC.Enh)) {
            thisHiC[j]<-GRanges(seqnames(thisHiC.TSS[j]),IRanges(
                min(start(thisHiC.TSS[j]),start(thisHiC.Enh[j])), ##select the minimum start position 
                max(end(thisHiC.TSS[j]),end(thisHiC.Enh[j]))))    ##select the maxmum end position 
            ##Make a CIGAR term for each
            thisHiC.cigar = paste0(width(thisHiC.Enh),"M",(width(thisHiC)-2*width(thisHiC.Enh)),"N",width(thisHiC.Enh),"M");
            #thisHiC.cigar<-paste0(width(thisHiC.Enh)%/%2,"M",(width(thisHiC)-2*width(thisHiC.Enh))+10000,"N",width(thisHiC.Enh)%/%2,"M");
            alTrack_HiC<-AlignmentsTrack(thisHiC,cigar=thisHiC.cigar,isPaired=FALSE,name=hic_name,type="sashimi", lwd.sashimiMax=3); # ,  col.anchors.fill ="#78AB41"
        }
        return(list(enhancer=thisHiC.Enh,tss=thisHiC.TSS,alTrack=alTrack_HiC))
    }
}

chromosome="chr11"  #Define the chromosome number
GABA<-Track_new(hicenhancer_GABA,hicpromoter_GABA,"GABA",chromosome)
Enh_GABA<-GABA$enhancer
tss_GABA<-GABA$tss
alTrack_GABA<-GABA$alTrack

GLU<-Track_new(hicenhancer_GLU,hicpromoter_GLU,"GLU",chromosome)
Enh_GLU<-GLU$enhancer
tss_GLU<-GLU$tss
alTrack_GLU<-GLU$alTrack

#window<-GRanges(seqnames=chromosome, 
#         IRanges(min(start(tss_pos),start(tss_neg),start(Enh_pos),start(Enh_neg))-200000, 
#                max(end(tss_pos),end(tss_neg),end(Enh_pos),end(Enh_neg))+200000))

#genomic region to plot

##Gene model
GeneModel<-BiomartGeneRegionTrack(genome="hg19", biomart=mart, chromosome=chromosome,
            start=120300000,end=121050000, showId=T, geneSymbols=T,
            transcriptAnnotation="symbol", name="Gene Model") 
generange<-GeneModel@range
#generange<-generange[generange$feature=='protein_coding']
GeneModel@range<-generange
GeneModel@dp@pars[1:25]<-"gray"
GeneModel@dp@pars$fill<-"gray"

htrack1<-HighlightTrack(trackList=GeneModel, start=unique(start(tss_GLU)), width=10000, chromosome=chromosome, fill="gray", alpha=0.5, col="#00000000")
#htrack3<-HighlightTrack(trackList=alTrack_GABA, start=unique(start(alTrack_GABA)), width=10000, chromosome=chromosome, fill="#339900", alpha=0.5, col="#00000000")
htrack4<-HighlightTrack(trackList=alTrack_GLU, start=unique(start(alTrack_GLU)), width=10000, chromosome=chromosome, fill="#993366", alpha=0.5, col="#00000000")

##track
tadtrack<-function(tadranges,sample_id,col_id){
    tadregions<-tadranges[seqnames(tadranges)==chromosome]
    tadtrack<-AnnotationTrack(tadregions, genome="hg19",name=sample_id, fill=col_id) 
    displayPars(tadtrack)<-list(alpha.title = 1, alpha = 0.5)
    return(tadtrack)
}

tad_track1<-tadtrack(tad_GABA_ran,"GABA TAD","#339900")
tad_track2<-tadtrack(tad_GLU_ran,"GLU TAD","#993366")

peak_track3<-tadtrack(GABAranges,"GABA peak","#339900")
peak_track4<-tadtrack(GLUranges,"GLU peak","#993366")

peak_track1<-DataTrack(range="/proj/hyejunglab/NeuN/ChIPseq/bam/SOX.27ac_200bp_CPM.bw",
            genome = "hg19",col="#339900",type="a",name="GABA peak", chromosome=chromosome)
peak_track2<-DataTrack(range="/proj/hyejunglab/NeuN/ChIPseq/bam/GLU.27ac_200bp_CPM.bw", 
            genome = "hg19",col="#993366",type="a",name="GLU peak",chromosome=chromosome)

gtrack<-GenomeAxisTrack(genome="hg19", chromosome=chromosome) # genome axis track
itrack<-itrackchr[[as.numeric(unlist(strsplit(chromosome, "chr"))[2])]] # chromosome figure

pdf(paste0("Loops_",genelist[1],".pdf"),width=7,height=7);
plotTracks(c(itrack,gtrack,htrack1,peak_track1,peak_track3,peak_track2,peak_track4,htrack4,tad_track1,tad_track2),shape="arrow",from=120300000,to=121050000,
    transcriptAnnotation="symbol",collapseTranscripts="meta",htrack4,add53=TRUE,showBandID=TRUE,cex.bands=0.5,stackHeight=0.5,
    background.title = "white",col.axis="black",col.title="black",cex.title=0.4,cex.axis=0.3,
    just.group="right",cex.main=0.8);  # 
dev.off()