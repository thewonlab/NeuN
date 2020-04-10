options(stringsAsFactors=F)
setwd("/proj/hyejunglab/crossdisorder/LDSC/sumstat/")
bashout = "/proj/hyejunglab/crossdisorder/LDSC/partherit/bash/"
directory = dir()
condition = unlist(strsplit(directory, split=".sumstats.gz"))
condition = c('ad.jansen','adhd2','asd2','bp2019','mdd.2019','pd2019','scz3','als','ms')

for(i in condition){
  system(paste0("echo '#!/bin/tcsh' > ",bashout, i,".tcsh"))
  system(paste0("echo '/proj/hyejunglab/program/ldsc/ldsc.py ",
                "--h2 /proj/hyejunglab/crossdisorder/LDSC/sumstat/", i,".sumstats.gz ",
                "--out /proj/hyejunglab/crossdisorder/LDSC/partherit/results/phase3/NeuN/", i, "_results ",
                "--frqfile-chr /proj/hyejunglab/program/ldsc/LDSC/1000G_Phase3_frq/1000G.EUR.QC. ",
                "--overlap-annot --ref-ld-chr /proj/hyejunglab/crossdisorder/LDSC/partherit/annot_NeuN/,/proj/hyejunglab/program/ldsc/LDSC/baselineLD_v1.1/baselineLD. ",
                "--w-ld-chr /proj/hyejunglab/program/ldsc/LDSC/weights_hm3_no_hla/weights.",
                "' >> ",bashout, i,".tcsh"))  
  system(paste0("sbatch -n 1 --mem=100g -o ", bashout, i, ".out ", bashout, i,".tcsh"))
}

condition =('smoking.GSCAN.2019','alcohol.GSCAN.2019','cannabis')
for(i in condition){
  system(paste0("echo '#!/bin/tcsh' > ",bashout, i,".tcsh"))
  system(paste0("echo '/proj/hyejunglab/program/ldsc/ldsc.py ",
                "--h2 /proj/hyejunglab/crossdisorder/LDSC/sumstat/", i,".sumstats.gz ",
                "--out /proj/hyejunglab/crossdisorder/LDSC/partherit/results/phase3/NeuN3/", i, "_results ",
                "--frqfile-chr /proj/hyejunglab/program/ldsc/LDSC/1000G_Phase3_frq/1000G.EUR.QC. ",
                "--overlap-annot --ref-ld-chr /proj/hyejunglab/crossdisorder/LDSC/partherit/annot_NeuN/,/proj/hyejunglab/program/ldsc/LDSC/baselineLD_v1.1/baselineLD. ",
                "--w-ld-chr /proj/hyejunglab/program/ldsc/LDSC/weights_hm3_no_hla/weights.",
                "' >> ",bashout, i,".tcsh"))  
  system(paste0("sbatch -n 1 --mem=100g -o ", bashout, i, ".out ", bashout, i,".tcsh"))
}