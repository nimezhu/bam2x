echo "    
# working dir  e.g. $HOME/scratch/job....
WORKING_DIR=.

# where you put your samples
# sample name (fastq file name without the suffix )
# if your fastq file is  SampleA.fastq  
# then  SAMPLE=SampleA
#       FASTQ_SUFFIX=fastq
#       PAIRED_END=0
# if you hava paired end data
#     your have two files ;  e.g.  SampleA_1.fastq SampleA_2.fastq
# then  SAMPLE=SampleA
#       FASTQ_SUFFIX=fastq
#       PAIRED_END=1
SAMPLE_DIR=.
SAMPLE=small_test     
FASTQ_SUFFIX=fastq
PAIRED_END=1          #1 or 0


# bwa index directory
# usally should put genome bwa_index file in nearline
#  the genome index should in the bwa_index directory
BWA_INDEX=$HOME/nearline/bwa_index
GENOME=mm10      

UNIQ_MAP=1   # 1 or 0  if it is 1 , it will report another bam which only have uniq map. 

# where you put the bam output file
OUTPUT_DIR=output

# where you submit the job to hpcc or just generate the scripts

SUBMIT=1   # 0 or 1 

MAIL=results.qsub@gmail.com


#=========================== USING THE DEFAULT  IF YOU WANT ===============================
# where you put samtools and bwa
BIN=/home/wangj2/bin                      
BWA=\$BIN/bwa
SAMTOOLS=\$BIN/samtools
# bwa aln options 
CPU= 100     # >24 will look for the available independent computer with max cpu.
BWA_ALN_OPTION=\"-Y\"
"


