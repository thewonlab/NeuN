###### super FIREs
cd /proj/hyejunglab/epigenetics/NeuNsort/Panos/Loops_ATAC
## analyze differential motifs with gimme maelstrom
echo -e "loc\tcluster"  > merged_peaks.txt

peak=(glia_loops_ATACpeak.txt neuron_loops_ATACpeak.txt)
sampleid=(NeuN- NeuN+)

for i in {0..1}
do
    input=${peak[$i]}
    cat $input | awk -v marker=${sampleid[$i]} 'NR>1 {print $1":"$2"-"$3"\t"marker}' >> merged_peaks.txt
done
conda create -n gimme python=3 gimmemotifs
conda activate gimme
genomepy install hg19 UCSC --annotation
gimme maelstrom merged_peaks.txt hg19 merged_peaks_motif.out

python
from gimmemotifs.maelstrom import MaelstromResult
import matplotlib.pyplot as plt
mr = MaelstromResult("merged_peaks_motif.out/")
fig = plt.figure(figsize=(10, 10))
mr.plot_heatmap(threshold=3)
plt.savefig("merged_peaks_motif.out/heatmap.pdf", bbox_inches = 'tight')