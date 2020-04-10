cd /proj/hyejunglab/NeuN/GWAS/
script=/nas/longleaf/home/hubenxia/project/GWAS/neuron_vs_glia/MAGAMA/SNP_position.py
dos2unix $script
snp=/proj/hyejunglab/disorder/Alzheimer/commonvar/Jansen_TableS8.csv #AD credible SNP list
neg_anno=/proj/hyejunglab/crossdisorder/annotation/NeuNneg_wointron.genes.annot
pos_anno=/proj/hyejunglab/crossdisorder/annotation/NeuNpos_wointron.genes.annot
python $script -i ${pos_anno} -r $snp -o AD_pos_SNP_position.txt
python $script -i ${neg_anno} -r $snp -o  AD_neg_SNP_position.txt

### run MAGMA

/proj/hyejunglab/crossdisorder/MAGMA/summarydat/

alzheimers.ctg.2018_MAGMA_Pval_sum.txt


module add magma
neg_anno=/proj/hyejunglab/crossdisorder/annotation/NeuNneg_wointron.genes.annot
pos_anno=/proj/hyejunglab/crossdisorder/annotation/NeuNpos_wointron.genes.annot
g1000_eur=/proj/hyejunglab/program/magma/ref_pop/g1000_eur/g1000_eur
cd /proj/hyejunglab/crossdisorder/MAGMA/summarydat/
arr=(alzheimers.ctg.2018_MAGMA_Pval_sum.txt)

for var in ${arr[@]}
do
    echo $var
    file_name=${var%%.*}
    magma --bfile $g1000_eur --pval $var use=rsid,P ncol=N --gene-annot $pos_anno --out /proj/hyejunglab/NeuN/GWAS/${file_name}_neuron
    magma --bfile $g1000_eur --pval $var use=rsid,P ncol=N --gene-annot $neg_anno --out /proj/hyejunglab/NeuN/GWAS/${file_name}_glia
done



cd /nas/longleaf/home/hubenxia/project/GWAS/neuron_vs_glia/MAGAMA
module add r/3.6.0
anno=/proj/hyejunglab/chr/geneAnno_allgenes.rda
outdir=/proj/hyejunglab/NeuN/GWAS/
Trajectory=/proj/hyejunglab/expression/kangetal/Kang_16874genes_1340samples_HJ.rda

input=("alzheimers")
disease=("AD")

Rscript --vanilla Developmental_expressionV2.R -a $anno -i ${input[0]} -d ${disease[0]} -t $Trajectory  -o /proj/hyejunglab/NeuN/GWAS/


### check gene expression in scRNA-seq
cd /nas/longleaf/home/hubenxia/project/GWAS/neuron_vs_glia/MAGAMA
module add r/3.6.0
scRNAseq=/proj/hyejunglab/singlecell/capstone4/lake2016quake2015/Capstone4_singlecellexpression_adultbrain.txt
Neuron='neuron'
Glia='glia'
outdir=/proj/hyejunglab/NeuN/GWAS/

neuron_ad=/proj/hyejunglab/NeuN/GWAS/AD_neuron_fdrthrgenes_adjusted.rda
glia_ad=/proj/hyejunglab/NeuN/GWAS/AD_glia_fdrthrgenes_adjusted.rda

neuron_gene=($neuron_ad)
glia_gene=($glia_ad)
disease=("AD")
Rscript --vanilla Cell_type_expression.R -n ${neuron_gene[0]} -g ${glia_gene[0i]} -N 'neuron' -G 'glia' -d ${disease[0]} -s $scRNAseq  -o /proj/hyejunglab/NeuN/GWAS/


neuron_ad=/proj/hyejunglab/NeuN/GWAS/AD_neuron_fdrthrgenes_adjusted.rda
glia_ad=/proj/hyejunglab/NeuN/GWAS/AD_glia_fdrthrgenes_adjusted.rda
neuron_gene=($neuron_ad)
glia_gene=($glia_ad)
disease=("AD")
Rscript --vanilla Cell_type_expressionV_AD.R -n ${neuron_gene} -g ${glia_gene} -N 'neuron' -G 'glia' -d ${disease} -s $scRNAseq  -o /proj/hyejunglab/NeuN/GWAS/

#plot venndiagram for AD Rrisk genes detected in neuron and glia
cd /nas/longleaf/home/hubenxia/project/GWAS/neuron_vs_glia/MAGAMA
module add r/3.6.0

neuron_ad=/proj/hyejunglab/NeuN/GWAS/AD_neuron_genename.txt
glia_ad=/proj/hyejunglab/NeuN/GWAS/AD_glia_genename.txt

neuron_gene=($neuron_ad)
glia_gene=($glia_ad)
disease=("AD")
Rscript --vanilla Venn.R -n ${neuron_gene[0]} -g ${glia_gene[0]} -o ${disease[0]}


# Gene ontology enrichment analysis of AD risk genes
cd /nas/longleaf/home/hubenxia/project/GWAS/neuron_vs_glia/MAGAMA
module add r/3.6.0

neuron_ad=/proj/hyejunglab/NeuN/GWAS/alzheimers_neuron.genes.out
glia_ad=/proj/hyejunglab/NeuN/GWAS/alzheimers_glia.genes.out

neuron_gene=($neuron_ad)
glia_gene=($glia_adz)
disease=("AD")
anno=/proj/hyejunglab/chr/geneAnno_allgenes.rda
Rscript --vanilla gProfileR.R -f ${neuron_gene[0]} -a $anno -o neuron_${disease[0i]}
Rscript --vanilla gProfileR.R -f ${glia_gene[0]} -a $anno -o glia_${disease[0]}


###extract protein-coding genes wit fdr<=0.05
cd /nas/longleaf/home/hubenxia/project/GWAS/neuron_vs_glia/MAGAMA
module add r/3.6.0
anno=/proj/hyejunglab/chr/geneAnno_allgenes.rda
outdir=/proj/hyejunglab/NeuN/GWAS/

input=("alzheimers")
disease=("AD")

Rscript --vanilla target_gene.R -a $anno -i ${input[0]} -d ${disease[0]}  -o /proj/hyejunglab/NeuN/GWAS/


# module analysis
cd /nas/longleaf/home/hubenxia/project/GWAS/neuron_vs_glia/MAGAMA
module add r/3.6.0
module=/proj/hyejunglab/NeuN/AD_H3K27ac/AD_disease/cels_208_mmc9.txt
anno=/proj/hyejunglab/chr/geneAnno_allgenes.rda
neuron_ad=/proj/hyejunglab/NeuN/GWAS/AD_neuron_genename.txt
glia_ad=/proj/hyejunglab/NeuN/GWAS/AD_glia_genename.txt
Rscript --vanilla AD_risk_gene_module.R -a $anno  -n ${neuron_ad} -g ${glia_ad} -m $module

# calculate gene expression for AD risk genes from scRNA-seq of AD.
module=/proj/hyejunglab/singlecell/AD/Mathys_Nature/
anno=/proj/hyejunglab/chr/geneAnno_allgenes.rda
neuron_ad=/proj/hyejunglab/NeuN/GWAS/AD_neuron_genename.txt
glia_ad=/proj/hyejunglab/NeuN/GWAS/AD_glia_genename.txt
Rscript --vanilla AD_risk_gene_scRNA_AD.R -a $anno  -n ${neuron_ad} -g ${glia_ad} -m $module -o 'AD_risk_gene'

