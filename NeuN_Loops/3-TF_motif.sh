## analyze differential motifs with gimme maelstrom
cd /proj/hyejunglab/epigenetics/NeuNsort/Panos/Loops_ATAC
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
