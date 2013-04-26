samtools mpileup -uf ../hg19.fna H3.sorted.bam | bcftools view -bvcg - > test_H3.var.raw.bcf
bcftools view test_H3.var.raw.bcf | vcfutils.pl varFilter -D100 > test_H3.var.flt.vcf
