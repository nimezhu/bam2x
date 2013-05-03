
echo "test average plot with two bams"
python plot_flanking_region.py -A -b /data/sysbio/MNase-seq_2013-3-18_mapping/Pn_E14_mm9_March13_rmdup.sort.bam /data/sysbio/MNase-seq_2013-3-18_mapping/Pn_KD_mm9_March13_rmdup.sort.bam -i CTCF_mm9_peak_5000.bed -f 150 -n E14 KD -o test_MNase_CTCF.png

echo "test heatmap plot with one bam and one interval"
python plot_flanking_region.py -H -b /data/sysbio/MNase-seq_2013-3-18_mapping/Pn_E14_mm9_March13_rmdup.sort.bam -i CTCF_mm9_peak_5000.bed -f 150 -o test_MNase_CTCF.png

echo "test heatmap plot with direction (two intervals)"
python plot_flanking_region.py -H -b /data/sysbio/MNase-seq_2013-3-18_mapping/Pn_E14_mm9_March13_rmdup.sort.bam -i test_mm9_5000_gene_1.bed test_mm9_5000_gene_2.bed -n mm9_gene-1 mm9_gene-2 -f 150 -o test_MNase_TSS.png -d -w 11

echo "#test average plot with direciotn (two intervals)"
python plot_flanking_region.py -A -b /data/sysbio/MNase-seq_2013-3-18_mapping/Pn_E14_mm9_March13_rmdup.sort.bam -i test_mm9_5000_gene_1.bed test_mm9_5000_gene_2.bed -n mm9_gene-1 mm9_gene-2 -f 150 -o test_MNase_TSS.png -d -w 11


