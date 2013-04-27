'''
Assume you start with a sorted.bam file of MNase-seq with index (input2)
'''
input2=MNase-E14_mm9_mapping.sorted.bam
echo "remove duplicates from sorted BAM file"
samtools rmdup -s $input2 MNase_E14_mm9_rmdup.sorted.bam
samtools index MNase_E14_mm9_rmdup.sorted.bam

input=MNase_E14_mm9_rmdup.sorted.bam
output=Pn_E14_rmdup_5e7.pdf

echo "Working on phasogram of $output with 3-pile"
python distogram_phasogram.py $input -b $input2 -n 50000000 -p -x 0 1500 -i 3 -o $output

'''
with the count file outputed from python scripts above "Phasogram_MNase_E14_rmdup_5e7pile-3_0~1500bp.txt", we can specify which region (for example (100,800) to plot by:
'''

python plot_histogram.py -p -x 100 800 -o new_phasogram.pdf
