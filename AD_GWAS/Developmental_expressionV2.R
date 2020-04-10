#!/usr/bin/env Rscript
library(optparse)
options(stringsAsFactors = FALSE)
library(biomaRt)
library(reshape)
library(WGCNA)
library(ggplot2)
library(gridExtra)
library(GenomicRanges)

option_list<-list(
    make_option(c("-a", "--anno"), type="character", default=NULL, 
            help="anno file name", metavar="character"),
    make_option(c("-i", "--input"), type="character", default=NULL, 
            help="input file name", metavar="character"),
    make_option(c("-d", "--disease"), type="character", default=NULL, 
            help="disease file name", metavar="character"),
    make_option(c("-t", "--trajectory"), type="character", default=NULL, 
            help="trajectory file name", metavar="character"),
    make_option(c("-o", "--outdir"), type="character", default=NULL, 
            help="outdir file name", metavar="character")
)

opt_parser<-OptionParser(option_list=option_list)
opt<-parse_args(opt_parser)

if (is.null(opt$anno)){
    print_help(opt_parser)
    stop("At least one argument must be supplied (anno file).n", call.=FALSE)
}
if (is.null(opt$input)){
    print_help(opt_parser)
    stop("At least one argument must be supplied (input file).n", call.=FALSE)
}
if (is.null(opt$disease)){
    print_help(opt_parser)
    stop("At least one argument must be supplied (disease file).n", call.=FALSE)
}
if (is.null(opt$trajectory)){
    print_help(opt_parser)
    stop("At least one argument must be supplied (trajectory name).n", call.=FALSE)
}
if (is.null(opt$outdir)){
    print_help(opt_parser)
    stop("At least one argument must be supplied (outdir file).n", call.=FALSE)
}

setwd(opt$outdir)
load(opt$anno) # Annotation for all genes 
pcgene = unique(geneAnno1[geneAnno1$gene_biotype=="protein_coding","ensembl_gene_id"])  # protein-coding genes
diseasename = c(opt$disease)

celltype = c("glia","neuron")
fdrthr= 0.05
mhcrange = GRanges(seqnames=6, IRanges(25000000,35000000)) # MHC region

getout<-function(filename,celltype,mhcrange,diseasename){
    diseasemat = read.table(paste0(filename,"_",celltype,".genes.out"), header=T) # reading MAGMA output files
    diseasepc = diseasemat[diseasemat$CHR!="X" & diseasemat$GENE %in% pcgene, ] # only get protein-coding genes
    # selecting out genes that are located in MHC region
    diseaseranges = GRanges(seqnames=diseasepc$CHR, IRanges(diseasepc$START, diseasepc$STOP), gene=diseasepc$GENE)
    olap = findOverlaps(diseaseranges, mhcrange)
    mhcgene = diseaseranges[queryHits(olap)]$gene
    diseasepc = diseasepc[!(diseasepc$GENE %in% mhcgene),] # remove MHC genes
    diseasepc$fdr = p.adjust(diseasepc$P, "BH") # calculate FDR from P-value
    selectgenes = diseasepc[diseasepc$fdr<fdrthr, "GENE"] # select ouf significantly associated genes
    diseasefdr = selectgenes
    genename = geneAnno1[geneAnno1$ensembl_gene_id %in% diseasefdr, "hgnc_symbol"]
    genename = genename[genename!=""]
    print(paste(filename,length(diseasefdr)))
    save(diseasefdr, genename,diseasepc, file=paste0(opt$outdir,diseasename,"_",celltype,"_fdrthrgenes_adjusted.rda"))
    write.table(diseasefdr, file=paste0(opt$outdir,diseasename,"_",celltype,'_geneid.txt'),quote=F,row.names=F,col.names=F)
    write.table(genename, file=paste0(opt$outdir,diseasename,"_",celltype,'_genename.txt'),quote=F,row.names=F,col.names=F)
    return(diseasefdr)
}

diseasefdr.glia<-getout(opt$input,celltype[1],mhcrange,diseasename)
diseasefdr.neuron<-getout(opt$input,celltype[2],mhcrange,diseasename)


## 1) Plot Expression Trajectory
load(opt$trajectory)
datExpr = scale(datExpr,center = T, scale=FALSE) + 1 # normalizing expression values 
datMeta$devstage = NA
datMeta[datMeta$age %in% c("5.7 PCW","6 PCW"),"devstage"] <- 1
datMeta[datMeta$age %in% c("8 PCW","9 PCW"),"devstage"] <- 2
datMeta[datMeta$age %in% c("12 PCW"),"devstage"] <- 3
datMeta[datMeta$age %in% c("13 PCW"),"devstage"] <- 4
datMeta[datMeta$age %in% c("16 PCW","17 PCW"),"devstage"] <- 5
datMeta[datMeta$age %in% c("19 PCW","21 PCW","22 PCW"),"devstage"] <- 6
datMeta[datMeta$age %in% c("25 PCW","35 PCW","37 PCW"),"devstage"] <- 7
datMeta[datMeta$age %in% c("4 PMonth"),"devstage"] <- 8
datMeta[datMeta$age %in% c("6 PMonth","10 PMonth"),"devstage"] <- 9
datMeta[datMeta$age %in% c("1 Y","2 Y","3 Y","4 Y"),"devstage"] <- 10
datMeta[datMeta$age %in% c("8 Y","11 Y"),"devstage"] <- 11
datMeta[datMeta$age %in% c("13 Y","15 Y","18 Y","19 Y"),"devstage"] <- 12
datMeta[datMeta$age %in% c("21 Y","22 Y","23 Y","27 Y","28 Y","30 Y","36 Y","37 Y","40 Y","42 Y","55 Y"),"devstage"] <- 13
datMeta[datMeta$age %in% c("64 Y","70 Y","82 Y"),"devstage"] <- 14


exprdat = vector(mode="list", length=length(celltype))

expr<-function(datExpr,target,datMeta,diseasename,filename){
    exprdat = vector(mode="list", length=length(diseasename))
    exprdat = apply(datExpr[match(target, rownames(datExpr)),],2,mean,na.rm=T)
    disname = diseasename
    dat = data.frame(Group=disname, Region=datMeta$Region, Age=datMeta$Age_unit, DevStage = datMeta$devstage, Unit=datMeta$Unit, Expr=exprdat, disname=diseasename)
    dat1 = dat
    dat1$Region = as.character(dat1$Region)
    dat1[dat1$Region %in% c("Frontal CTX","Temporal CTX","Parietal CTX","Occipital CTX"), "Region"] <- "CTX"
    dat1$Region = factor(dat1$Region)
    regionstoplot = c("CTX")
    SUD = diseasename
    disgroup = c(SUD)
    names(disgroup) = c(filename)
    pdf(file=paste0(opt$outdir,diseasename,"_developmental.regional.expression.fdr_",filename,"_union_adjusted.pdf"), width=6, height=5)
    for(i in 1:length(regionstoplot)){
        dat2 = dat1[dat1$Region==regionstoplot[i],]
        dat2 = dat2[dat2$Group %in% disgroup,]
        print(head(dat2))
        p.brainspanTime = ggplot(dat2,aes(x=DevStage, y= Expr, fill=Group, color=Group)) + xlim(1,15) + ylab("Normalized expression") + 
            geom_smooth(span=1) + xlab("Developmental Stages") + scale_x_continuous(breaks=seq(1,14,1))+
            ggtitle(paste0("Brain Developmental Expression Trajectory:",regionstoplot[i])) 
        print(p.brainspanTime)
    }
    dev.off()

    dat.sud = dat1[dat1$Group %in% diseasename,]
    print(head(dat.sud))
    p.brainExp.sud = 
        ggplot(dat.sud,aes(x=Unit, y=Expr, fill=Unit)) + 
        ylab("Normalized expression") +  geom_boxplot(outlier.size = NA)  +
        ggtitle("Brain Expression") + xlab("") +
        scale_alpha_manual(values=c(0.2, 1))

    pdf(paste0(opt$outdir,diseasename,"_prenatalVSpostnatal_union_adjusted_",filename,".pdf"), height=5, width=8)
    #print(p.brainExp.sud)
    dev.off()
    dat1$Index = paste(dat1$Group, dat$Unit, sep="_")
    dat1$Index = unlist(lapply(strsplit(as.character(dat1$Index), split=" "),'[[',1))
    dat1$Index = factor(dat1$Index)
    a<-dat1[dat1$Index=='AD_Prenatal','Expr']
    b<-dat1[dat1$Index=='AD_Postnatal','Expr']
    print(head(as.numeric(a)))
    print(head(as.numeric(b)))
    print(wilcox.test(as.numeric(a),as.numeric(b))$p.value)
    fit <- aov(Expr ~ Index, dat1)
    summary(fit)
    pval = TukeyHSD(fit)$Index
    pvaltoreport = paste0(disname, "_Prenatal (Post-Conception Week)-",disname,"_Postnatal (Years)")
    pvaltoreport = pval[rownames(pval) %in% pvaltoreport,]
    fitsum = c()
    compareset1 = paste(disname,"Prenatal",sep="_")
    compareset2 = paste(disname,"Postnatal",sep="_")
    compareset = dat1[dat1$Index %in% c(compareset1, compareset2),]
    fit = lm(Expr ~ Index, compareset)
    fitsum = rbind(fitsum,summary(fit)$coefficients[2,])
    rownames(fitsum) = disname
    fitsum = data.frame(fitsum)
    colnames(fitsum) = c("Estimate", "SE","t.value","P")
    return(fitsum)
}
glia_statistics<-expr(datExpr,diseasefdr.glia,datMeta,diseasename,celltype[1])
neuron_statistics<-expr(datExpr,diseasefdr.neuron,datMeta,diseasename,celltype[2])
print(neuron_statistics)
print(glia_statistics)