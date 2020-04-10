library(reshape)
library(ggplot2)

dat<-data.frame(
regulation=c("T-M1","T-M5","T-M16","T-M9","T-M4","T-M7","T-M8","T-M18","T-M10","T-M3","T-M14"),
M1=c(2,0,0,0,0,0,0,0,0,0,0),
M5=c(0,2,0,0,0,0,0,0,0,0,0),
M16=c(0,0,2,0,0,0,0,0,0,0,0),
M9=c(0,0,0,2,0,0,0,0,0,0,0),
M4=c(0,0,0,0,2,0,0,0,0,0,0),
M7=c(0,0,0,0,0,1,0,0,0,0,0),
M8=c(0,0,0,0,0,0,1,0,0,0,0),
M18=c(0,0,0,0,0,0,0,1,0,0,0),
M10=c(0,0,0,0,0,0,0,0,1,0,0),
M3=c(0,0,0,0,0,0,0,0,0,1,0),
M14=c(0,0,0,0,0,0,0,0,0,0,1)
)

level<-c("T-M1","T-M5","T-M16","T-M9","T-M4","T-M7","T-M8","T-M18","T-M10","T-M3","T-M14")
melt.data<-melt(dat, id.vars="regulation", variable_name="flag")
pdf("module_regulation.pdf")
ggplot(data=melt.data,aes(x=factor(regulation,levels=level),y=flag,fill=factor(value)))+geom_tile(color = "black")+
    scale_fill_manual(values=c("0"="white", "1"="red","2"="blue"))+
    labs(y="Module", x="") +theme_minimal() +
    theme_bw() + theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank())
dev.off()

dat<-data.frame(
regulation=c("T-M1","T-M5","T-M16","T-M9","T-M4","T-M7","T-M8","T-M18","T-M10","T-M3","T-M14"),
Neuron=c(2,2,2,2,2,0,0,0,0,0,0),
Astrocytes=c(0,0,0,0,0,1,1,0,1,0,1),
Microglia=c(0,0,0,0,0,0,0,1,0,1,0),
Oligo=c(0,0,0,0,0,0,0,0,0,0,0),
Endothelial=c(0,0,0,0,0,0,0,0,0,0,0)
)

level<-c("T-M1","T-M5","T-M16","T-M9","T-M4","T-M7","T-M8","T-M18","T-M10","T-M3","T-M14")
melt.data<-melt(dat, id.vars="regulation", variable_name="flag")
pdf("module_regulation2.pdf")
ggplot(data=melt.data,aes(y=factor(regulation,levels=level),x=flag,fill=factor(value)))+geom_tile(color = "black")+
    scale_fill_manual(values=c("0"="white", "1"="red","2"="blue"))+
    labs(y="Module", x="") +theme_minimal() +
    theme_bw() + theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank())
dev.off()
