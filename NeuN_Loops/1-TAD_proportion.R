library(ggplot2)
options(stringsAsFactors=F)
setwd("/proj/hyejunglab/NeuN/Loop")

system(paste0("intersectBed -a NeuNneg_TSS_interacting_region_fixed.bed -b /proj/hyejunglab/NeuN/NeuNpos/chrmat/40kb/HMMdomain.bed -wa -wb > NeuNneg_TAD_evolloci.bed"))
system(paste0("intersectBed -a NeuNpos_TSS_interacting_region_fixed.bed -b /proj/hyejunglab/NeuN/NeuNneg/chrmat/40kb/HMMdomain.bed -wa -wb > NeuNpos_TAD_evolloci.bed"))

TADN<-read.table("NeuNneg_TAD_evolloci.bed")
TADP<-read.table("NeuNpos_TAD_evolloci.bed")

TADin<-c()
for(i in 1:dim(TADN)[1]){
  if(TADN[i,4]>TADN[i,7] & TADN[i,4]<TADN[i,8]){
    TADin<-c(TADin,TRUE)
  }else{
    TADin<-c(TADin,FALSE)
  }
}
TADN<-cbind(TADN, TADin)

TADin<-c()
for(i in 1:dim(TADP)[1]){
  if(TADP[i,4]>TADP[i,7] & TADP[i,4]<TADP[i,8]){
    TADin<-c(TADin,TRUE)
  }else{
    TADin<-c(TADin,FALSE)
  }
}
TADP<-cbind(TADP, TADin)
colnames(TADN)<-colnames(TADP)<-c("chr", "start", "end", "intstart", "intend", "tadchr", "tadstart", "tadend", 'NA','NA','NA','NA','NA','NA',"TADin")

pieN<-table(TADN$TADin)
names(pieN)<-c("not in the same TAD", "in the same TAD")
dfN = data.frame(tad = names(pieN), counts=as.vector(pieN))
pieN=ggplot(data=dfN, aes(x="",y=counts,fill=tad))

pieP = table(TADP$TADin); names(pieP) = c("not in the same TAD", "in the same TAD")
dfP = data.frame(tad = names(pieP), counts=as.vector(pieP))
pieP=ggplot(data=dfP, aes(x="",y=counts,fill=tad))

blank_theme <- theme_minimal()+
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.border = element_blank(),
    panel.grid=element_blank(),
    axis.ticks = element_blank(),
    plot.title=element_text(size=14, face="bold")
  )

pdf("TADin_TSS.pdf", width=10, height=5)
par(mfrow=c(1,2))
pieN + geom_bar(width=1, stat="identity")+coord_polar(theta="y") + scale_fill_hue(h=c(180,270)) + 
  theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank()) + blank_theme +
  geom_text(aes(label = counts)) 
pieP + geom_bar(width=1, stat="identity")+coord_polar(theta="y") + scale_fill_hue(h=c(180,270)) + 
  theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank()) + blank_theme +
  geom_text(aes(label = counts)) 
dev.off()

piedat = rbind(dfN, dfP)
piedat$celltype=c("NeuN-", "NeuN-", "NeuN+", "NeuN+")
pdf("TADin_barplot_TSS.pdf", width=5, height=5)
barpl = ggplot(data=piedat, aes(x=celltype,y=counts, fill=as.factor(tad)))
barpl + geom_bar(position="fill", stat="identity") + theme_bw() + 
  theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(), axis.text.x=element_text(angle=90,hjust=1)) 
dev.off()
