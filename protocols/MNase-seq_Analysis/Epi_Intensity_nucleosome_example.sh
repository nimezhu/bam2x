# change ChIP-seq fragment length as 300

folder=~/ChIPseq_map_data/

nucleosome_file=/home/yu68/MNase-seq/Danpos/whole_MNase/result/pooled/Pn_E14_whole_mm9.sort.Fnor.ajClonal.smooth.peaks.xls

python Epi_Intensity_nucleosome.py -l 300 -N $nucleosome_file -b ${folder}mouse_H3K4me3_d0.bed ${folder}mouse_H3K4me2_d0.bed ${folder}mouse_H3K4me1_d0.bed ${folder}mouse_H3K27me3_d0.bed -n H3K4me3 H3K4me2 H3K4me1 H3K27me3 -o Epi_Intensity_nucleosome-1-300bp.txt

python Epi_Intensity_nucleosome.py -l 300 -N $nucleosome_file -b ${folder}mouse_H3K36me3_d0.bed ${folder}mouse_H3K9me3.bed ${folder}mouse_H3K27ac_d0.bed ${folder}mouse_H2AZ_d0.bed -n H3K36me3 H3K9me3 H3K27ac H2AZ -o Epi_Intensity_nucleosome-2-300bp.txt


