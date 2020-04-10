title: "NeuN_Alzheimers_CellTypeExpression"

options(stringsAsFactors=F)
library(ggplot2)
library(reshape2)
library(corrplot)
library(dplyr)
library(plyr)

setwd("/nas/longleaf/home/hubenxia/project/Alzheimer/neuron_vs_glia")
load("pos_hypo.genelist.rda")
write.table(file="pos_hypo.genelist.txt",as.data.frame(pos_hypo.genelist),quote=F,row.names=F)
load("neg_hypo.genelist.rda")
write.table(file="neg_hypo.genelist.txt",as.data.frame(neg_hypo.genelist),quote=F,row.names=F)
load("pos_hyper.genelist.rda")
write.table(file="pos_hyper.genelist.txt",as.data.frame(pos_hyper.genelist),quote=F,row.names=F)
load("neg_hyper.genelist.rda")
write.table(file="neg_hyper.genelist.txt",as.data.frame(neg_hyper.genelist),quote=F,row.names=F)


load("geneAnno_allgenes.rda")


#NeuNpos_hypo and NeuNneg_hyper only
cellexp<-read.table("/nas/longleaf/home/hubenxia/project/Loops/neuron_vs_glia/Capstone4_singlecellexpression_adultbrain.txt", header=T)
rownames(cellexp)<-cellexp[,1]
cellexp<-cellexp[,-1]
datExpr<-scale(cellexp,center=T, scale=F) 

theme_set(theme_bw())

summarySE<-function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
    library(plyr)
  # New version of length which can handle NA's: if na.rm==T, don't count them
    length2<-function(x, na.rm=FALSE) {
        if (na.rm){
            sum(!is.na(x))
        }else{
            length(x)
        }
    }
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
    datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
    #Rename the "mean" column    
    datac <- rename(datac, c("mean" = measurevar))
    datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
    #Confidence interval multiplier for standard error
    #Calculate t-statistic for confidence interval: 
    #e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
    ciMult <- qt(conf.interval/2 + .5, datac$N-1)
    datac$ci <- datac$se * ciMult
    return(datac)
}


target<-function(data1,data2,name1,name2){
    targetensg<-list(data1,data2)
    targetname<-c("NeuNpos_hypo", "NeuNneg_hyper")
    targetgene<-vector(length=length(targetensg), mode="list")
    for(i in 1:length(targetgene)){
        genelist<-targetensg[[i]]
        genelist<-geneAnno1[match(genelist,geneAnno1$hgnc_symbol),"hgnc_symbol"]
        genelist<-unique(genelist)
        genelist<-genelist[genelist!=""]
        genelist<-genelist[!is.na(genelist)]
        targetgene[[i]] <- genelist
    }

    exprdat<-vector(mode="list", length=length(targetgene))

    for(i in 1:length(targetgene)){
        exprdat[[i]]<-apply(datExpr[match(targetgene[[i]], rownames(datExpr)),],2,mean,na.rm=T)
    }

    dat<-c()
    for(i in 1:length(targetgene)){
        datframe = data.frame(Group=targetname[i], cell=names(exprdat[[i]]), Expr=exprdat[[i]])
        dat = rbind(dat, datframe)
    }

    dat$cell<-gsub("Fetal.quiescent","Fetal-quiescent",dat$cell)
    dat$cell<-gsub("Fetal.replicating","Fetal-replicating",dat$cell)
    dat$celltype<- unlist(lapply(strsplit(dat$cell, split="[.]"),'[[',1))
    dat$celltype<-factor(dat$celltype, levels=c("Neurons","Astrocytes","Microglia","Endothelial","Oligodendrocytes",                                        "OPC","Fetal-quiescent","Fetal-replicating",
                                "Ex1","Ex2","Ex3","Ex4","Ex5","Ex6","Ex7","Ex8",
                                "In1","In2","In3","In4","In5","In6","In7","In8"))

    sedat<-summarySE(dat, measurevar="Expr", groupvars=c("celltype", "Group"))
    scale_this <- function(x) as.vector(scale(x))
    sedat$subtype<-ifelse(sedat$celltype %in% c(paste0("Ex",1:8), paste0("In",1:8)),"subtype","main")
    library(dplyr)
    scaledat = sedat %>% 
        group_by(Group, subtype) %>% 
        mutate(value=scale_this(Expr))

    return(list(dat,sedat))
}


exp_ac<-function(sedate,type_class,group1,group2,filename){
    sedate1<-sedate[sedate$subtype==type_class & sedate$Group==group1,]
    sedate2<-sedate[sedate$subtype==type_class & sedate$Group==group2,]
    if(nrow(sedate1) & nrow(sedate2)){
        sedate1$value<- scale(sedate1$Expr)
        sedate2$value<- scale(sedate2$Expr)
        scalemate<-rbind(as.numeric(sedate1$value), as.numeric(sedate2$value))
        rownames(scalemate)<-c("NeuNpos_hypo", "NeuNneg_hyper")
        colnames(scalemate)<-sedate1$celltype
        scalemate<-scalemate
    }else if (!nrow(sedate1) & nrow(sedate2)){
        sedate2$value<- scale(sedate2$Expr)
        scalemate<-t(as.data.frame(as.numeric(sedate2$value)))
        rownames(scalemate)<-"NeuNneg_hyper"
        colnames(scalemate)<-sedate2$celltype
        scalemate<-scalemate
    }else if (nrow(sedate1) & !nrow(sedate2)){
        sedate1$value<- scale(sedate1$Expr)
        scalemate<-t(as.data.frame(as.numeric(sedate1$value)))
        rownames(scalemate)<-"NeuNpos_hyper"
        colnames(scalemate)<-sedate1$celltype
        scalemate<-scalemate
    }
    print(scalemate)
    pdf(paste0(filename,"_singlecell_exp_ind_",type_class,".pdf"), width=10, height=5)
    col2plot = colorRampPalette(c("#053061","#2166AC","#4393C3","#92C5DE","#D1E5F0","#FFFFFF","#FDDBC7","#F4A582","#D6604D","#B2182B","#67001F"))
    corrplot(scalemate, is.corr=FALSE, method="color", tl.col="black", col=col2plot(100)) # 
    dev.off()
}


#NeuNpos_hypo and NeuNneg_hyper only
sedat_both<-list()
exp_ac<-function(sedate,type_class,group1,group2,filename){
    sedate1<-sedate[sedate$subtype==type_class & sedate$Group==group1,]
    sedate2<-sedate[sedate$subtype==type_class & sedate$Group==group2,]
    if(nrow(sedate1) & nrow(sedate2)){
        sedate1$value<- scale(sedate1$Expr)
        sedate2$value<- scale(sedate2$Expr)
        scalemate<-rbind(as.numeric(sedate1$value), as.numeric(sedate2$value))
        rownames(scalemate)<-c("NeuNpos_hypo", "NeuNneg_hyper")
        colnames(scalemate)<-sedate1$celltype
    }
    print(scalemate)
    pdf(paste0(filename,"_singlecell_exp_ind_",type_class,".pdf"), width=10, height=5)
    col2plot = colorRampPalette(c("#053061","#2166AC","#4393C3","#92C5DE","#D1E5F0","#FFFFFF","#FDDBC7","#F4A582","#D6604D","#B2182B","#67001F"))
    corrplot(scalemate, is.corr=FALSE, method="color", tl.col="black", col=col2plot(100), cl.lim=c(-2.0,2.0)) # 
    dev.off()
}
sedat_both<-target(pos_hypo.genelist,neg_hyper.genelist,"NeuNpos_hypo", "NeuNneg_hyper")[[2]]
exp_ac(sedat_both,"main","NeuNpos_hypo", "NeuNneg_hyper","poshypo_neghyper")
exp_ac(sedat_both,"subtype","NeuNpos_hypo", "NeuNneg_hyper","poshypo_neghyper")

#Hypo-Cell Type Expression
sedat_Hypo<-target(pos_hypo.genelist, neg_hypo.genelist,"NeuNpos_hypo", "NeuNneg_hypo")[[2]]
exp_ac(sedat_Hypo,"main","NeuNpos_hypo", "NeuNneg_hypo","Hypoacetylated")
exp_ac(sedat_Hypo,"subtype","NeuNpos_hypo", "NeuNneg_hypo","Hypoacetylated")

#Hyper-Cell Type Expression
sedat_Hyper<-target(pos_hyper.genelist, neg_hyper.genelist,"NeuNpos_hyper", "NeuNneg_hyper")[[2]]
exp_ac(sedat_Hyper,"main","NeuNpos_hyper", "NeuNneg_hyper","Hyperacetylated")
exp_ac(sedat_Hyper,"subtype","NeuNpos_hyper", "NeuNneg_hyper","Hyperacetylated")
