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
