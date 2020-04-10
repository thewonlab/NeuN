## GWAS analysis pipeline
cd /nas/longleaf/home/hubenxia/project/Loops/GLU_GABA
module add r/3.6.0
GABA=/proj/hyejunglab/NeuN/ChIPseq/diffbind/GABA_diff_peaks.txt # GABA differential H3K27ac peaks 
Glu=/proj/hyejunglab/NeuN/ChIPseq/diffbind/Glu_diff_peaks.txt  # Glu differential H3K27ac peaks 
snp=/proj/hyejunglab/disorder/Schizophrenia_pgc/imputedsnp/allsnp.bed
Rscript --vanilla file_for_magma.R -s $snp -p $GABA -o 'GABA_SCZ'
Rscript --vanilla file_for_magma.R -s $snp -p $Glu -o 'Glu_SCZ'

refer=/proj/hyejunglab/crossdisorder/annotation/NeuNpos_wointron.genes.annot
dos2unix long2wide.py
python long2wide.py -i GABA_SCZ_snp_long.txt -r $refer -o GABA_SCZ_.genes.annot
python long2wide.py -i Glu_SCZ_snp_long.txt -r $refer -o Glu_SCZ_.genes.annot

script=/nas/longleaf/home/hubenxia/project/Loops/GLU_GABA/SNP_position.py
dos2unix $script
snp=/proj/hyejunglab/disorder/Schizophrenia_pgc/imputedsnp/allsnp.bed
python $script -i GABA_SCZ_.genes.annot -r $snp -o GABA_SNP_position.txt
python $script -i Glu_SCZ_.genes.annot -r $snp -o  Glu_SNP_position.txt

module add r/3.6.0
R
Glu<-read.table("Glu_SCZ_snp_long.txt",sep="\t")
GABA<-read.table("GABA_SCZ_snp_long.txt",sep="\t")
Glu_peak<-nrow(unique(Glu[,1:3]))   #50901
GABA_peak<-nrow(unique(GABA[,1:3])) #37661

Glu_snp<-length(unique(Glu[,4]))   #368025
GABA_snp<-length(unique(GABA[,4])) #306154


### run MAGMA

module add magma
GABA_anno=/proj/hyejunglab/NeuN/Glu_GABA/MAGMA/GABA_SCZ_.genes.annot
Glu_anno=/proj/hyejunglab/NeuN/Glu_GABA/MAGMA/Glu_SCZ_.genes.annot
g1000_eur=/proj/hyejunglab/program/magma/ref_pop/g1000_eur/g1000_eur
cd /proj/hyejunglab/crossdisorder/MAGMA/summarydat/
arr=(scz.clozuk.2017_MAGMA_Pval_sum.txt BD.PGC.2019_MAGMA_Pval_sum.txt)
for var in ${arr[@]}
do
    echo $var
    file_name=${var%%.*}
    magma --bfile $g1000_eur --pval $var use=rsid,P ncol=N --gene-annot $GABA_anno --out /proj/hyejunglab/NeuN/Glu_GABA/MAGMA/${file_name}_GABA
    magma --bfile $g1000_eur --pval $var use=rsid,P ncol=N --gene-annot $Glu_anno --out /proj/hyejunglab/NeuN/Glu_GABA/MAGMA/${file_name}_Glu
done


## convert geneid to gene symbol on the basis of FDR less than 0.05
script=/nas/longleaf/home/hubenxia/project/Loops/GLU_GABA/0.05/target_gene.R
module add r/3.6.0
anno=/proj/hyejunglab/chr/geneAnno_allgenes.rda
cd /proj/hyejunglab/NeuN/Glu_GABA/MAGMA/

input=("BD" "scz")
disease=("BD" "SCZ")

for i in {0..1}
do
    Rscript --vanilla $script -a $anno -i ${input[$i]} -d ${disease[$i]}  -o /nas/longleaf/home/hubenxia/project/Loops/GLU_GABA/0.05/
done

## cell type expression profile

cd /nas/longleaf/home/hubenxia/project/Loops/GLU_GABA/0.05/
module add r/3.6.0
scRNAseq=/proj/hyejunglab/singlecell/capstone4/lake2016quake2015/Capstone4_singlecellexpression_adultbrain.txt
outdir=/nas/longleaf/home/hubenxia/project/Loops/GLU_GABA/0.05/
GABA_scz=/nas/longleaf/home/hubenxia/project/Loops/GLU_GABA/0.05/SCZ_GABA_fdrthrgenes_adjusted.rda
Glu_scz=/nas/longleaf/home/hubenxia/project/Loops/GLU_GABA/0.05/SCZ_Glu_fdrthrgenes_adjusted.rda
GABA_bd=/nas/longleaf/home/hubenxia/project/Loops/GLU_GABA/0.05/BD_GABA_fdrthrgenes_adjusted.rda
Glu_bd=/nas/longleaf/home/hubenxia/project/Loops/GLU_GABA/0.05/BD_Glu_fdrthrgenes_adjusted.rda
GABA_gene=($GABA_bd $GABA_scz)
Glu_gene=($Glu_bd $Glu_scz)
disease=("BD" "SCZ")
for i in {0..1}
do
    Rscript --vanilla Cell_type_expression.R -S ${GABA_gene[$i]} -B ${Glu_gene[$i]} -N 'GABA' -G 'Glu' -d ${disease[$i]} -s $scRNAseq  -o /nas/longleaf/home/hubenxia/project/Loops/GLU_GABA/0.05/
done