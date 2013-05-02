#!/bin/sh

inputfile=$1

# you can add friendly code to check the extension and use zcat if it's *.fastq.gz or *.fq.gz...

[[ "$inputfile" =~ ".*\.(fq\.gz$|fastq\.gz$)" ]] && zcat $inputfile | head -n400 | awk '{if(NR%4==0) printf("%s",$0);}' |  od -A n -t u1 | awk 'BEGIN{min=100;max=0;}{for(i=1;i<=NF;i++) {if($i>max) max=$i; if($i<min) min=$i;}}END{if(max<=74 && min<59) print "Phred+33"; else if(max>73 && min>=64) print "Phred+64"; else if(min>=59 && min<64 && max>73) print "Solexa+64"; else print "Unknown score encoding!";}'


[[ "$inputfile" =~ ".*\.(fq$|fastq$)" ]] && head -n400 $inputfile | awk '{if(NR%4==0) printf("%s",$0);}' |  od -A n -t u1 | awk 'BEGIN{min=100;max=0;}{for(i=1;i<=NF;i++) {if($i>max) max=$i; if($i<min) min=$i;}}END{if(max<=74 && min<59) print "Phred+33"; else if(max>73 && min>=64) print "Phred+64"; else if(min>=59 && min<64 && max>73) print "Solexa+64"; else print "Unknown score encoding!";}'
