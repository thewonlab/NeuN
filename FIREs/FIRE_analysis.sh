##### FIRE analysis for neuronal and glial Hi-C
cd /proj/hyejunglab/NeuN/FIRE
#### call FIREcaller to analyze FIREs and superFIREs for neuronal and glial Hi-C
### run 1-FIRE.R to obtain FIREs and superFIREs.

#### analyze differential FIREs
Rscript --vanilla differentialFIREs.R  -f FIRE_ANALYSIS_40KB.txt  -N 'neuron' -G 'glia'

Rscript --vanilla differential_FIREs_plot.R  -f FIRE_ANALYSIS_40KB.txt  -o 'neuron_glia'

#### differential FIREs overlap with differential H3K27ac peaks
neuron=/nas/longleaf/home/hubenxia/project/FIRE/neuronvsglia/neuron_H3K27ac_differential_peak.txt
glia=/nas/longleaf/home/hubenxia/project/FIRE/neuronvsglia/non-neuron_H3K27ac_differential_peak.txt
Rscript --vanilla differential_FIRE_overlap_diffH3K27ac.R  -f FIRE_ANALYSIS_40KB.txt  -n $neuron -g $glia -N 'neuron' -G 'glia'

#########################we did it for paper
#### FIREs overlap with promoters
cd /proj/hyejunglab/NeuN/FIRE
anno=/proj/hyejunglab/NeuN/FIRE/geneAnno_allgenes.rda
promoter=/proj/hyejunglab/NeuN/FIRE/gencode19_promoter.bed
Rscript --vanilla differentialFIREs_overlap_promoter.R -p $promoter -f neuron_diffFIREs.bed -g $anno -o neuron
Rscript --vanilla differentialFIREs_overlap_promoter.R -p $promoter -f glia_diffFIREs.bed -g $anno -o glia

#### GO analysis with gProfileR

Rscript --vanilla gProfile2.R  -f neuron_diffFIRE_promoter_gene.txt  -o 'neurondiffFIRE_gene'
Rscript --vanilla gProfile2.R  -f glia_diffFIRE_promoter_gene.txt  -o 'gliadiffFIRE_gene'

scRNA=/proj/hyejunglab/NeuN/FIRE/Capstone4_singlecellexpression_adultbrain.txt
neuron=/proj/hyejunglab/NeuN/FIRE/neuron_diffFIRE_promoter_gene.txt
glia=/proj/hyejunglab/NeuN/FIRE/glia_diffFIRE_promoter_gene.txt
Rscript --vanilla FIRE_NeuNposVSneg.R  -f $scRNA  -n $neuron -g $glia -o 'no'
