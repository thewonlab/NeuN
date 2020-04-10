cd /proj/hyejunglab/GLUGABA
out=/proj/hyejunglab/NeuN/ChIPseq
for line in `ls *.fastq.gz`
do
   fastqc  $line -o $out
done

cd /proj/hyejunglab/GLUGABA

arr=(H.276 H.286 H.344 H.372 H.395 H.406 H.412 H.427 H.444)
sample=(GLU.27ac GLU.input SOX.27ac SOX.input)
ref=/proj/seq/data/HG19_UCSC/Sequence/Bowtie2Index/genome
out=/proj/hyejunglab/NeuN/ChIPseq/bam
for var in ${arr[@]}
do
    for varid in ${sample[@]}
    do
        left=${var}.${varid}.R1.fastq.gz
        right=${var}.${varid}.R2.fastq.gz
        file_name=${var}.${varid}
        bowtie2 --very-sensitive -p 16 -x  $ref -1 $left  -2 $right -S $out/${file_name}.sam
        samtools view -bS $out/${file_name}.sam  > $out/${file_name}.bam
        rm $out/${file_name}.sam
        samtools sort $out/${file_name}.bam -o $out/${file_name}_sort.bam
        rm $out/${file_name}.bam
        picard MarkDuplicates I=$out/${file_name}_sort.bam O=$out/${file_name}_clean.bam M=$out/dups.txt REMOVE_DUPLICATES=true
        rm $out/${file_name}_sort.bam
        samtools sort $out/${file_name}_clean.bam -o $out/${file_name}.bam
        rm $out/${file_name}_clean.bam
        samtools index $out/${file_name}.bam
    done
done

cd /proj/hyejunglab/NeuN/ChIPseq/bam

arr=(H.276 H.286 H.344 H.372 H.395 H.406 H.412 H.427 H.444)
cell=(GLU SOX)
#IP=(27ac input)
num=(0 1)
out=/proj/hyejunglab/NeuN/ChIPseq/macs2
for var in ${arr[@]}
do
    for i in ${num[@]}
    do
        IP=${var}.${cell[$i]}.27ac.bam
        input=${var}.${cell[$i]}.input.bam
        macs2 callpeak -t $IP -c $input --outdir $out -g hs -n ${var}.${cell[$i]} -B -f BAM --broad --broad-cutoff 0.00001 
    done
done

cd /proj/hyejunglab/NeuN/ChIPseq/bam
cell=(GLU SOX)
IP=(27ac input)
num=(0 1)
for i in ${num[@]}
do
    for j in ${num[@]}
    do
        samtools merge ${cell[$i]}.${IP[$j]}.bam *.${cell[$i]}.${IP[$j]}.bam
        samtools index ${cell[$i]}.${IP[$j]}.bam
    done
done

cell=(GLU SOX)
IP=(27ac input)
num=(0 1)
for i in ${num[@]}
do
    bamCompare -b1 ${cell[i]}.${IP[0]}.bam -b2 ${cell[i]}.${IP[1]}.bam -o ${cell[i]}.${IP[0]}_200bp_CPM.bw --operation ratio -bs 200 -p 8 --normalizeUsing CPM  --scaleFactorsMethod None
done


### analyze differential peaks with diffbind
run diffbind.R

### generate differential peaks in Glu and GABA cells
cd /proj/hyejunglab/NeuN/ChIPseq/diffbind
cat H3K27ac3_differential_peaks.txt | awk 'OFS="\t"{if($9>0) print $0}' > Glu_diff_peaks.txt  #56037
cat H3K27ac3_differential_peaks.txt | awk 'OFS="\t"{if($9<0) print $0}' > GABA_diff_peaks.txt  #41951


cd /nas/longleaf/home/hubenxia/project/ChIP-seq/GLU_GABA
module add r/3.6.0
dir=/proj/hyejunglab/NeuN/ChIPseq/macs2
Rscript --vanilla  mergePeaks.R -d $dir -p GLU -o  /proj/hyejunglab/NeuN/ChIPseq/macs2/Glu_peak.bed
Rscript --vanilla  mergePeaks.R -d $dir -p GABA -o  /proj/hyejunglab/NeuN/ChIPseq/macs2/GABA_peak.bed
