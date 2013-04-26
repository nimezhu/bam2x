#for i in chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY chrM
for i in   chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chrX chrY chrM
do
echo "samtools mpileup -uf ../hg19.fna ../Mapping/$i.v2.sorted.bam | bcftools view -bvcg - > test_$i.var.raw.bcf"
echo "bcftools view test_$i.var.raw.bcf | vcfutils.pl varFilter -D100000 > test_$i.vcf"
echo "vcf2var.pl test_$i.vcf > test_$i.variants"
echo "samtools view ../Mapping/$i.v2.sorted.bam | extract_hairs --sam stdin --variants test_$i.variants > test_$i.fragments"
done 
