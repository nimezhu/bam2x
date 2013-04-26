 for i in H3K27ac H3K36me3 H3K4me1 H3K4me3 H3K9me3
  do
 #     for j in chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY

#     do
# ./xHapBlockBinomialTest.py ./HAP/test_$j.haplotype ../Mapping/${i}_H209_hg19_mapping_v2.sorted.bam ./NCI-H209.copynumber.tab> ./HapBlockBinomialTest/${i}.HapBlockBinomialTest.tab.$j
#  grep Simple ./HapBlockBinomialTest/$i.HapBlockBinomialTest.tab.$j > ./HapBlockBinomialTest/Simple/$i.HapBlockBinomialTest.tab.$j.simple
#  
#
# done
 cat ./HapBlockBinomialTest/Simple/$i.*.simple > ./HapBlockBinomialTest/Raw/$i.HapBlockBinomialTest.raw.tab
 done
