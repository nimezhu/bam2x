#!/bin/sh

#  STEP0. reading the config file
#  Using config file with the same filename with different suffix
#  e.g. script.sh using script.config 

# config file example:
# ======== START =====================
# WORKING_DIR=.
# BWA_INDEX=$HOME/nearline/bwa_index
# BIN=/home/wangj2/bin                      
# SAMPLE_DIR=.
# OUTPUT_DIR=output
# SAMPLE=small_test
# CPU=8
# GENOME=mm10
# SAMTOOLS=$BIN/samtools
# BWA=$BIN/bwa
# BWA_ALN_OPTION="-Y"
# MAIL=yourname@gmail.com
# FASTQ_SUFFIX=fastq           # (fastq or fq)
# PAIRED_END=1                 # (1 or 0)
# =========== END ========================

function usage()
{
    echo "
    USAGE: $0 CONFIG_FILE 
           or
           $0   ( using default config file ${0%.*}.config )


    "
}


if [ $# -ge 1 ] && [ -f $1 ]
    then
        CONFIG_FILE=$1
    else
        CONFIG_FILE=${0%.*}.config
fi


if [ ! -f $CONFIG_FILE ]
    then
        echo "ERROR: Can't find the config file $CONFIG_FILE"
        usage
        exit

fi

. $CONFIG_FILE

echo $CPU
if [ $CPU -gt 24 ]
    then
    CPU=`qcpu.sh -m`
fi

echo $CPU

#=================================================================================================
# STEP1. mapping two paired end fastq file (*_[12].fq or *_[12].fastq) using bwa



STEP1=${WORKING_DIR}/bwa_aln.${SAMPLE}_${GENOME}.qsub.step1
echo "
#$ -V
#$ -cwd
#$ -pe single $CPU
#$ -o $HOME/sge_jobs_output/sge_job.\$JOB_ID.out -j y        
#$ -S /bin/bash
#$ -m abe
#$ -M $MAIL

" > $STEP1




if [ $PAIRED_END -eq 1 ]
    then
echo "
$BWA aln "$BWA_ALN_OPTION" -t $CPU $BWA_INDEX/${GENOME}.fa ${SAMPLE_DIR}/${SAMPLE}_1.$FASTQ_SUFFIX -f ${WORKING_DIR}/${SAMPLE}_1.aln.sai 
$BWA aln "$BWA_ALN_OPTION" -t $CPU $BWA_INDEX/${GENOME}.fa ${SAMPLE_DIR}/${SAMPLE}_2.$FASTQ_SUFFIX -f ${WORKING_DIR}/${SAMPLE}_2.aln.sai 
mv $HOME/sge_jobs_output/sge_job.\$JOB_ID.out $OUTPUT_DIR/bwa_aln.${SAMPLE}_${GENOME}.qsub.step1.log
"   >> ${WORKING_DIR}/bwa_aln.${SAMPLE}_${GENOME}.qsub.step1
    else
echo "
$BWA aln "$BWA_ALN_OPTION" -t $CPU $BWA_INDEX/${GENOME}.fa ${SAMPLE_DIR}/${SAMPLE}.$FASTQ_SUFFIX -f ${WORKING_DIR}/${SAMPLE}.aln.sai

if [ ! -d $OUTPUT_DIR ]
    then
    mkdir $OUTPUT_DIR
fi
if [ ! -d $OUTPUT_DIR/log ]
    then
    mkdir $OUTPUT_DIR/log
fi


mv $HOME/sge_jobs_output/sge_job.\$JOB_ID.out $OUTPUT_DIR/log/bwa_aln.${SAMPLE}_${GENOME}.qsub.step1.log
"   >> $STEP1
 
fi




#==========================================================================================
#  STEP2. one thread , using fewer CPUs.
#  . merge paired end alignments file into sam
#  . transfer sam file into bam format
#  . remove the sam file and the two alignment file *.sai
#  . sort bam file
#  . remove the unsorted bam file
#  . index the sorted bam file
#  . output the sorted bam file and its index into output directory

STEP2=${WORKING_DIR}/bwa_aln.${SAMPLE}_${GENOME}.qsub.step2
echo "
#$ -V
#$ -cwd
#$ -pe single 4
#$ -o $HOME/sge_jobs_output/sge_job.\$JOB_ID.out -j y        
#$ -S /bin/bash
#$ -m abe
#$ -M $MAIL
" > $STEP2


if [ $PAIRED_END -eq 1 ]
    then
echo "
$BWA sampe -P $BWA_INDEX/${GENOME}.fa ${WORKING_DIR}/${SAMPLE}_1.aln.sai ${WORKING_DIR}/${SAMPLE}_2.aln.sai ${SAMPLE_DIR}/${SAMPLE}_1.fastq ${SAMPLE_DIR}/${SAMPLE}_2.fastq -f ${WORKING_DIR}/${SAMPLE}_${GENOME}.sam
$SAMTOOLS view -S -b -o ${WORKING_DIR}/${SAMPLE}_${GENOME}.bam ${WORKING_DIR}/${SAMPLE}_${GENOME}.sam
rm ${WORKING_DIR}/${SAMPLE}_1.aln.sai 
rm ${WORKING_DIR}/${SAMPLE}_2.aln.sai
" >> $STEP2

    else

echo "

$BWA samse $BWA_INDEX/${GENOME}.fa ${WORKING_DIR}/${SAMPLE}.aln.sai ${SAMPLE_DIR}/${SAMPLE}.fastq  -f ${WORKING_DIR}/${SAMPLE}_${GENOME}.sam
$SAMTOOLS view -S -b -o ${WORKING_DIR}/${SAMPLE}_${GENOME}.bam ${WORKING_DIR}/${SAMPLE}_${GENOME}.sam
rm ${WORKING_DIR}/${SAMPLE}.aln.sai 
" >>$STEP2
fi


echo "
rm ${WORKING_DIR}/${SAMPLE}_${GENOME}.sam
$SAMTOOLS sort ${WORKING_DIR}/${SAMPLE}_${GENOME}.bam ${WORKING_DIR}/${SAMPLE}_${GENOME}.sorted
rm ${WORKING_DIR}/${SAMPLE}_${GENOME}.bam
$SAMTOOLS index ${WORKING_DIR}/${SAMPLE}_${GENOME}.sorted.bam
" >> $STEP2



### UNIQ MAPPING
if [ $UNIQ_MAP -eq 1 ]
    then
echo "
$SAMTOOLS view -h ${WORKING_DIR}/${SAMPLE}_${GENOME}.sorted.bam | egrep \"^@|XT:A:U\"  | $SAMTOOLS view -Sb - > ${WORKING_DIR}/${SAMPLE}_${GENOME}.uniqmap.sorted.bam 
$SAMTOOLS index ${WORKING_DIR}/${SAMPLE}_${GENOME}.uniqmap.sorted.bam 
" >> $STEP2
fi



echo "
if [ ! -d $OUTPUT_DIR ]
    then
    mkdir $OUTPUT_DIR
fi
if [ ! -d $OUTPUT_DIR/log ]
    then
    mkdir $OUTPUT_DIR/log
fi


if [ \"$OUTPUT_DIR\"!=\"$WORKING_DIR\" ]
    then
      mv ${WORKING_DIR}/${SAMPLE}_${GENOME}.* $OUTPUT_DIR/
fi

mv $HOME/sge_jobs_output/sge_job.\$JOB_ID.out $OUTPUT_DIR/log/bwa_aln.${SAMPLE}_${GENOME}.qsub.step2.log
" >> $STEP2



#===========================================================================
# STEP 3 PUT THE SCRIPTS into OUTPUT_DIR/log

if [ ! -d $OUTPUT_DIR ]
    then
    mkdir $OUTPUT_DIR
fi


if [ ! -d $OUTPUT_DIR/log ]
    then
    mkdir $OUTPUT_DIR/log
fi





if [ "$OUTPUT_DIR"!="$WORKING_DIR" ]
    then
        cp ${WORKING_DIR}/bwa_aln.${SAMPLE}_${GENOME}.qsub.step1  ${OUTPUT_DIR}/log/
        cp ${WORKING_DIR}/bwa_aln.${SAMPLE}_${GENOME}.qsub.step2  ${OUTPUT_DIR}/log/
fi


#=========================================================================
# STEP 4 SUBMIT THE qsub in dependency.


#if [ $SUBMIT -eq 1 ]
#    then
#        job1=`qsub $STEP1`
#        job2=`qsub -W depend=afterok:$job1 $STEP2`
#    else
#        echo "
#        job1=\`qsub $STEP1\`
#        job2=\`qsub -W depend=afterok:\$job1id $STEP2\`
#        "  > ${WORKING_DIR}/bwa_aln.${SAMPLE}_${GENOME}.qsub.exe
#        chmod 755 ${WORKING_DIR}/bwa_aln.${SAMPLE}_${GENOME}.qsub.exe
#fi

if [ $SUBMIT -eq 1 ]
    then
        job1=`qsub $STEP1`
        job1_id=`echo $job1 | cut -f 3 -d ' '`
        job2=`qsub -hold_jid $job1_id $STEP2`
    else
        echo "
        job1=\`qsub $STEP1\`
        job1_id=\`echo \$job1 | cut -f 3 -d ' '\`
        job2=\`qsub -hold_jid \$job1_id $STEP2\`

        "  > ${WORKING_DIR}/bwa_aln.${SAMPLE}_${GENOME}.qsub.exe
        chmod 755 ${WORKING_DIR}/bwa_aln.${SAMPLE}_${GENOME}.qsub.exe
fi



