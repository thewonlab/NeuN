#Installing EWCE:
install.packages("devtools")
library(devtools)
install.packages("ggdendro")
install_github("nathanskene/ewce")

#You can then load the package:
library(EWCE)
library(ggplot2)
library(cowplot)
library(limma)
library(readxl)

#Loading datasets
cellexp = read.table("/nas/longleaf/home/hubenxia/project/FIRE/neuronvsglia/Capstone4_singlecellexpression_adultbrain.txt", header=T)
rownames(cellexp) = cellexp[,1]
cellexp = cellexp[,-1]


cellexp = fix.bad.hgnc.symbols(cellexp,dropNonHGNC =TRUE)
annotation<-read.table("/nas/longleaf/home/hubenxia/project/FIRE/neuronvsglia/scRNA-seq_annotation.txt",sep="\t",header=T)
level1class=as.vector(annotation$level_1)
level2class=as.vector(annotation$level_2)

exp_CortexOnly_DROPPED = drop.uninformative.genes(
exp=as.matrix(cellexp),level2annot = level2class)

annotLevels = list(level1class=level1class,level2class=level2class)

sctd = generate.celltype.data(
exp=exp_CortexOnly_DROPPED,annotLevels=annotLevels,
groupName="sctd")

sctd<-load("/nas/longleaf/home/hubenxia/project/FIRE/neuronvsglia/FIRE/CellTypeData_sctd.rda")
human.bg = read.table("/nas/longleaf/home/hubenxia/project/FIRE/neuronvsglia/HGNC.txt",header=T,stringsAsFactors=F)[,1]
#human.bg = unique(human.bg)
#Application to genetic data
#Preparing gene lists
setwd("/proj/hyejunglab/NeuN/FIRE")
posgene = read.table("neuron_diffFIRE_promoter_gene.txt", header=F,sep="\t")  ##no extend
neggene = read.table("glia_diffFIRE_promoter_gene.txt", header=F,sep="\t")  ##no extend


posgene_genelist<-as.vector(posgene[,1])
neggene_genelist<-as.vector(neggene[,1])

reps=15000

human.hits = posgene_genelist
human.bg = unique(c(human.hits,human.bg))

pos_results =bootstrap.enrichment.test(sct_data=ctd,hits=human.hits, bg=human.bg,reps=reps,annotLevel=1,geneSizeControl=TRUE,genelistSpecies="human",sctSpecies="human")

human.hits = neggene_genelist
human.bg = unique(c(human.hits,human.bg))
neg_results = bootstrap.enrichment.test(sct_data=ctd,hits=human.hits,
                    bg=human.bg,reps=reps,annotLevel=1,geneSizeControl=TRUE,genelistSpecies="human",sctSpecies="human")
cell_name=c("Ex","In","Neurons","Astrocytes","Microglia",
"Oligodendrocytes","OPC","Endothelial",
"Fetal.quiescent","Fetal.replicating")

pos_results$results$CellType <- factor(pos_results$results$CellType,levels= cell_name,labels= cell_name)
neg_results$results$CellType <- factor(neg_results$results$CellType,levels= cell_name,labels= cell_name)

full_res2 = data.frame(pos_results$results,list="NeuNpos")
scnd_res2 = data.frame(neg_results$results,list="NeuNneg")
merged_results = rbind(full_res2,scnd_res2)


##############FIRE-extend
#pdf("/nas/longleaf/home/hubenxia/project/FIRE/neuronvsglia/FIRE/significant_results_cell_type-level1_control_FIRE_GREAT.pdf")
#ewce.plot(total_res=merged_results,mtc_method="BH")
#dev.off()

##############FIRE-no
pdf("/nas/longleaf/home/hubenxia/project/FIRE/neuronvsglia/FIRE/significant_results_cell_type-level1_control_FIRE_GREAT_no.pdf")
ewce.plot(total_res=merged_results,mtc_method="BH")
dev.off()