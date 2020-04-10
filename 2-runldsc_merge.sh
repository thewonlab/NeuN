#!/bin/bash
#SBATCH -p general
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --time=5-00:00:00
#SBATCH --mem=500g
#SBATCH --mail-type=end
#SBATCH --mail-user=hubenxia@gmail.com
#SBATCH -c 16

cat_dir="/proj/hyejunglab/crossdisorder/LDSC/partherit/annot_Glu_GABA" ## This is the directory for the annotation file
cd $cat_dir
for chr in `seq 1 22`
do
    python /proj/hyejunglab/program/ldsc/ldsc.py --l2 --bfile  /proj/hyejunglab/program/ldsc/LDSC/1000G_EUR_Phase3_plink/1000G.EUR.QC.$chr  --ld-wind-cm 1 --annot $cat_dir/$chr.annot.gz --out $cat_dir/$chr --print-snps    /proj/hyejunglab/program/ldsc/LDSC/list.txt
done



