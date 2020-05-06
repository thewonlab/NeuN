##### FIRE analysis for neuronal and glial Hi-C
cd /proj/hyejunglab/NeuN/FIRE
#### call FIREcaller to analyze FIREs and superFIREs for neuronal and glial Hi-C
### run 1-FIRE.R to obtain FIREs and superFIREs.

#### analyze differential FIREs
Rscript --vanilla differentialFIREs.R  -f FIRE_ANALYSIS_40KB.txt  -N 'neuron' -G 'glia'

#### differential FIREs overlap with differential H3K27ac peaks
neuron=/nas/longleaf/home/hubenxia/project/FIRE/neuronvsglia/neuron_H3K27ac_differential_peak.txt
glia=/nas/longleaf/home/hubenxia/project/FIRE/neuronvsglia/non-neuron_H3K27ac_differential_peak.txt
Rscript --vanilla differential_FIRE_overlap_diffH3K27ac.R  -f FIRE_ANALYSIS_40KB.txt  -n $neuron -g $glia -N 'neuron' -G 'glia'

#### cell expression of genes assigned to differential FIREs
scRNA=/proj/hyejunglab/NeuN/FIRE/Capstone4_singlecellexpression_adultbrain.txt
neuron=/proj/hyejunglab/NeuN/FIRE/neuron_diffFIRE_promoter_gene.txt
glia=/proj/hyejunglab/NeuN/FIRE/glia_diffFIRE_promoter_gene.txt
Rscript --vanilla FIRE_NeuNposVSneg.R  -f $scRNA  -n $neuron -g $glia -o 'no'
