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
options(stringsAsFactors = F)

setwd("/nas/longleaf/home/hubenxia/project/Alzheimer/neuron_vs_glia")
neunpos_hypo<-read.table("hypo_neuron_overlap.txt",header=T,sep="\t")
neunneg_hyper<-read.table("hyper_glia_overlap.txt",sep="\t",header=T)

enhrange<-function(input_file){
    enhancer_range<-GRanges(seqnames=input_file[,1],IRanges(input_file[,2],input_file[,3]),tss=input_file[,6],gene=input_file[,8])
    return(enhancer_range)
}

prorange<-function(input_file){
    promoter_range<-GRanges(seqnames=input_file[,1],IRanges(input_file[,6],input_file[,7]),tss=input_file[,6],gene=input_file[,8])
    return(promoter_range)
}

hicenhancer_pos<-enhrange(neunpos_hypo)
hicpromoter_pos<-prorange(neunpos_hypo)
hicenhancer_neg<-enhrange(neunneg_hyper)
hicpromoter_neg<-prorange(neunneg_hyper)

##Generate other required files
tad_pos<-read.table("/proj/hyejunglab/NeuN/NeuNpos/chrmat/HMMdomain.bed")
tad_neg<-read.table("/proj/hyejunglab/NeuN/NeuNneg/chrmat/HMMdomain.bed")
tad_range<-function(tad){
    tadranges<-GRanges(seqnames=tad[,1], IRanges(tad[,2], tad[,3]))
    return(tadranges)
}
tad_pos_ran<-tad_range(tad_pos)
tad_neg_ran<-tad_range(tad_neg)

############all differential peaks
neuron<-read.table("neuron_H3K27ac_differential_peak.txt",sep="\t",header=T,quote="")[,c(2:4)]
glia<-read.table("non-neuron_H3K27ac_differential_peak.txt",sep="\t",header=T,quote="")[,c(2:4)]
neuronranges<-tad_range(neuron)
gliaranges<-tad_range(glia)
hypo <- read.csv("AD_hypoacetylated.csv",header=T)
hyper <- read.csv("AD_hyperacetylated.csv",header=T)
granges <- function(input){
    input$CHR[input$CHR == 23]<- "X"
    input$CHR[input$CHR == 24] <- "Y"
    result <- GRanges(seqnames=paste0('chr',input$CHR), IRanges(input$start, input$end))
    return(result)
}
hyporanges<-granges(hypo)
hyperranges<-granges(hyper)

genelist<-c("CACNG3")
chromosome=unique(neunpos_hypo[neunpos_hypo[,8]==genelist,1])
peak_start=unique(neunpos_hypo[neunpos_hypo[,8]==genelist,9])
peak_end=unique(neunpos_hypo[neunpos_hypo[,8]==genelist,10])

Control <- read.table("/proj/hyejunglab/NeuN/NeuN_Alzheimer_acetylation/H3K27ac_and_HiC/conunionS.tab",sep="\t",skip=1)
AD <- read.table("/proj/hyejunglab/NeuN/NeuN_Alzheimer_acetylation/H3K27ac_and_HiC/ADunionS.tab",sep="\t",skip=1)
Control[is.na(Control)] <- 0
AD[is.na(AD)] <- 0
Control <- GRanges(seqnames=Control[,1], ranges=IRanges(start=Control[,2],end=Control[,3]),value=Control[,4])
AD <- GRanges(seqnames=AD[,1], ranges=IRanges(start=AD[,2],end=AD[,3]),value=AD[,4])


peak_track3<-DataTrack(range=Control,
            genome = "hg19",col="#339900",type="a",name="NeuNpos peak", chromosome=chromosome)
peak_track1<-DataTrack(range="/proj/hyejunglab/NeuN/Loop/ABC/H3K27ac_neuron_200b_RPKM.bw",ylim=c(0, 20),
            genome = "hg19",col="#339900",type="a",name="NeuNpos peak", chromosome=chromosome)
peak_track4<-DataTrack(range=AD,
            genome = "hg19",col="#993366",type="a",name="NeuNneg peak",chromosome=chromosome)
peak_track2<-DataTrack(range="/proj/hyejunglab/NeuN/Loop/ABC/H3K27ac_glia_200b_RPKM.bw", ylim=c(0, 20),
            genome = "hg19",col="#993366",type="a",name="NeuNneg peak",chromosome=chromosome)

load("/proj/hyejunglab/chr/geneAnno_allgenes.rda")
load("/proj/hyejunglab/chr/itrackchromosome.rda")

mart <- useMart(biomart="ENSEMBL_MART_ENSEMBL",dataset="hsapiens_gene_ensembl",host="feb2014.archive.ensembl.org")

Track_new<-function(enh_range,pro_range,hic_name,genelist,chr_name){
    chrnum<-chr_name # chromosome number
    thisHiC.Enh<-enh_range
    thisHiC.TSS<-pro_range
    thisHiC.Enh<-thisHiC.Enh[thisHiC.Enh$gene==genelist]  ##we need to ensure thatthisHiC.Enh is not NULL
    thisHiC.TSS<-thisHiC.TSS[thisHiC.TSS$gene==genelist]  ##we need to ensure thatthisHiC.TSS is not NULL
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
            print(thisHiC[j])
            thisHiC.cigar = paste0(width(thisHiC.Enh),"M",(width(thisHiC)-2*width(thisHiC.Enh)),"N",width(thisHiC.Enh),"M");
            #thisHiC.cigar<-paste0(width(thisHiC.Enh)%/%2,"M",(width(thisHiC)-2*width(thisHiC.Enh))+10000,"N",width(thisHiC.Enh)%/%2,"M");
            print(thisHiC.cigar)
            alTrack_HiC<-AlignmentsTrack(thisHiC,cigar=thisHiC.cigar,isPaired=FALSE,name=hic_name,type="sashimi", lwd.sashimiMax=3); # ,  col.anchors.fill ="#78AB41"
        }
        return(list(enhancer=thisHiC.Enh,tss=thisHiC.TSS,alTrack=alTrack_HiC))
    }
}

Pos<-Track_new(hicenhancer_pos,hicpromoter_pos,"NeuNpos",genelist,chromosome)
Enh_pos<-Pos$enhancer
tss_pos<-Pos$tss
alTrack_pos<-Pos$alTrack

Neg<-Track_new(hicenhancer_neg,hicpromoter_neg,"NeuNneg",genelist,chromosome)
Enh_neg<-Neg$enhancer
tss_neg<-Neg$tss
alTrack_neg<-Neg$alTrack

#window<-GRanges(seqnames=chromosome, 
#         IRanges(min(start(tss_pos),min(start(tss_neg),start(Enh_pos),start(Enh_neg))-200000, 
#                max(end(tss_pos),max(end(tss_neg),end(Enh_pos),end(Enh_neg))+200000))
window<-GRanges(seqnames=chromosome, 
         IRanges(min(start(tss_pos),start(Enh_pos))-100000, 
                max(end(tss_pos),end(Enh_pos))+150000))
##Gene model
GeneModel<-BiomartGeneRegionTrack(genome="hg19", biomart=mart, chromosome=chromosome,
            start=start(window), end=end(window), showId=T, geneSymbols=T,
            transcriptAnnotation="symbol", name="Gene Model") 
generange<-GeneModel@range
generange<-generange[generange$feature=='protein_coding',]
GeneModel@range<-generange
GeneModel@dp@pars[1:25]<-"gray"
GeneModel@dp@pars$fill<-"gray"

htrack1<-HighlightTrack(trackList=GeneModel, start=unique(start(tss_pos)), width=10000, chromosome=chromosome, fill="gray", alpha=0.5, col="#00000000")
htrack3<-HighlightTrack(trackList=alTrack_pos, start=unique(start(alTrack_pos)), width=10000, chromosome=chromosome, fill="#339900", alpha=0.5, col="#00000000")
#htrack4<-HighlightTrack(trackList=alTrack_neg, start=unique(start(alTrack_neg)), width=10000, chromosome=chromosome, fill="#993366", alpha=0.5, col="#00000000")

tadtrack<-function(tadranges,sample_id,col_id){
    tadregions<-tadranges[seqnames(tadranges)==chromosome]
    tadtrack<-AnnotationTrack(tadregions, genome="hg19",name=sample_id, fill=col_id) 
    displayPars(tadtrack)<-list(alpha.title = 1, alpha = 0.5)
    return(tadtrack)
}

tad_track1<-tadtrack(tad_pos_ran,"NeuNpos TAD","#339900")
tad_track2<-tadtrack(tad_neg_ran,"NeuNneg TAD","#993366")

peak_track3<-tadtrack(neuronranges,"NeuNpos peak","#339900")
peak_track4<-tadtrack(gliaranges,"NeuNneg peak","#993366")
peak_track7<-tadtrack(hyporanges,"hypo peak","#339900")
peak_track8<-tadtrack(hyperranges,"hyper peak","#993366")

peak_track5<-DataTrack(range=Control,
            genome = "hg19",col="#339900",type="a",name="NeuNpos peak", chromosome=chromosome)
peak_track1<-DataTrack(range="/proj/hyejunglab/NeuN/Loop/ABC/H3K27ac_neuron_200b_RPKM.bw",ylim=c(0, 20),
            genome = "hg19",col="#339900",type="a",name="NeuNpos peak", chromosome=chromosome)
peak_track6<-DataTrack(range=AD,
            genome = "hg19",col="#993366",type="a",name="NeuNneg peak",chromosome=chromosome)
peak_track2<-DataTrack(range="/proj/hyejunglab/NeuN/Loop/ABC/H3K27ac_glia_200b_RPKM.bw", ylim=c(0, 20),
            genome = "hg19",col="#993366",type="a",name="NeuNneg peak",chromosome=chromosome)


gtrack<-GenomeAxisTrack(genome="hg19", chromosome=chromosome) # genome axis track
itrack<-itrackchr[[as.numeric(unlist(strsplit(chromosome, "chr"))[2])]] # chromosome figure
ht<-HighlightTrack(trackList = c(peak_track1,peak_track3,peak_track5,peak_track7,peak_track2,peak_track4,peak_track6,peak_track8) ,start=peak_start, end=peak_end, chromosome = chromosome)

pdf(paste0("Loops_",genelist,".pdf"),width=7,height=7)
plotTracks(c(itrack,gtrack,htrack1,ht,htrack3,tad_track1,tad_track2),from=start(window),to=end(window),
    transcriptAnnotation="symbol",collapseTranscripts="meta",add53=TRUE,showBandID=TRUE,cex.bands=0.5,stackHeight=0.5,
    background.title = "white",col.axis="black",col.title="black",cex.title=0.4,cex.axis=0.3,
    just.group="right",cex.main=0.8);  # collapseTranscripts="meta"
dev.off()

