library(Vennerable)
fire<-read.table("FIRE_ANALYSIS_40KB.txt", header=T)
fire$counts<-fire$NeuNpos_indicator + fire$NeuNneg_indicator
pos<-fire[fire$NeuNpos_indicator==1,][,1:3]
neg<-fire[fire$NeuNneg_indicator==1,][,1:3]
pos_vector<-paste(pos[,1],pos[,2],pos[,3],sep="-")
neg_vector<-paste(neg[,1],neg[,2],neg[,3],sep="-")
FIRE<-list(pos_FIRE=pos_vector,neg_FIRE=neg_vector)
Vdemo<-Venn(FIRE)
pdf("FIRE_overlap.pdf")
plot(Vdemo, doWeights = TRUE, type = "circles")
dev.off()