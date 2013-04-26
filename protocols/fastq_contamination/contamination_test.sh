#!/bin/bash

if [ $# -lt 1 ]
    then
    echo "USAGE: test_contamination.sh <fastq_file>"
    exit 0
fi

# modify the following parameters as needed
Sample_Size=5000000
Fastq_to_fasta_PATH=/home/zocean/Software/fastx_toolkit-0.0.13/src/fastq_to_fasta/fastq_to_fasta
BLAST_PATH=/home/zocean/Software/ncbi-blast-2.2.26+/bin/blastn

# sample fastq file
./fastqsampler.py $1 $Sample_Size >${1%.fastq}_sample_$Sample_Size.fastq

# convert fastq to fasta
$Fastq_to_fasta_PATH -n -v -i ${1%.fastq}_sample_$Sample_Size.fastq  -o ${1%.fastq}_sample_$Sample_Size.fasta

## download Salmo Salar genome sequence from ncbi or you can download na nr database
## wget http://www.ncbi.nlm.nih.gov/Traces/wgs/?download=AGKD01.fasta.gz
## reference assembly papge
## http://www.ncbi.nlm.nih.gov/genome/assembly/313068/

## build Salmo salar genome blast detabase
## /home/zocean/Software/ncbi-blast-2.2.26+/bin/makeblastdb -in 2011-10-25_Salmo_salar_WGS_AGKD01.fasta -parse_seqids -hash_index -out Salmo_salar_blast_database_AGKD01 -taxid 8030 -title Salmo_salar_WGS -input_type fasta -dbtype nucl

# do blast
$BLAST_PATH -query ${1%.fastq}_sample_$Sample_Size.fasta -db Salmo_salar_blast_database_AGKD01 -out blast.result -max_target_seqs 5 -num_threads 4 -outfmt "7 qacc sacc length evalue score pident" -lcase_masking


sed '/# BLAST/d' blast_result.txt | sed '/# Data/d' |sed '/# Field/d' | sed 's/^#\s//' | awk '{if(/Query/) {Num++;printf("%d\t%s",Num,$0);} else if(/DBRHHJN1/){printf("%d\t.\t.\t%s\n",Num,$0);} else {printf("%d\t%s\n",Num,$0);}}' >blast_parse.txt

awk 'BEGIN{num=0;sign=0}{if(sign==1) {print $0;sign=0;} if($1!=num) {num=$1;sign=1};}' blast_parse.txt | cut -f 5 | sed '/^$/d' | sed 's/^/gb|/;s/$/|/'> top1_ID.txt

# retrieve species ID
./retrieve_gb_by_acc.pl top1_ID.txt >top1_ID_retrieve_species.txt
cut -f 2 top1_ID_retrieve_species.txt | sort -k 1 | uniq -c | sed 's/\s\{1,\}//;s/\s/\t/' | sort -k 1 -n -r >pie_chart.txt
