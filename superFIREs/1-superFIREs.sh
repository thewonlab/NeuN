#!/usr/bin/env Rscript
library("FIREcaller")
library(MASS)
library(preprocessCore)
library(data.table)
prefix.list <- c('NeuNpos','NeuNneg')
gb<-'hg19'
map_file<-'F_GC_M_Hind3_40Kb_el.hg19.txt.gz'
rm_mhc <- TRUE
FIREcaller(prefix.list, gb, map_file, rm_mhc)
