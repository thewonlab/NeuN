library(Vennerable)
pos<-read.table("super_FIRE_call_NeuNpos.txt", header=T,sep="\t")[,1:3]
neg<-read.table("super_FIRE_call_NeuNneg.txt", header=T,sep="\t")[,1:3]
pos_vector<-paste(pos[,1],pos[,2],pos[,3],sep="-")
neg_vector<-paste(neg[,1],neg[,2],neg[,3],sep="-")
super_FIRE<-list(pos_superFIRE=pos_vector,neg_superFIRE=neg_vector)
Vdemo<-Venn(super_FIRE)
pdf("super_FIRE_overlap.pdf")
plot(Vdemo, doWeights = TRUE, type = "circles")
dev.off()

pos_gene<-read.table("NeuNposgene_superFIRE_H3K27ac.txt", header=T)  #gene list
neg_gene<-read.table("NeuNneggene_superFIRE_H3K27ac.txt", header=T)  #gene list
super_FIRE_gene<-list(pos_superFIRE=pos_gene,neg_superFIRE=neg_gene)
Vdemo<-Venn(super_FIRE_gene)
pdf("super_FIRE_gene_overlap.pdf")
plot(Vdemo, doWeights = TRUE, type = "circles")
dev.off()