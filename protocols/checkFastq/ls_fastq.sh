#/usr/bin/env sh

if [ $# -ge 1 ]
    then
      dir=$1
    else
      dir="./"
fi

for i0 in $dir/*.fastq
    do
       # dir=${i0%/*.fastq}
        i=${i0##*/}
        if [ "$i" == "*.fastq" ]
            then
              echo "no fastq in this directory"
              break
        fi
        j2=${i%_1.fastq}
        mate2="${j2}_2.fastq"
        j1=${i%_2.fastq}
        mate1="${j1}_1.fastq"
       if [ ! -f "$dir/$mate1" ]
       then
        if [ -f "$dir/$mate2" ]
            then
               echo "paired end $i $mate2"
            else
              echo "single end $i"
        fi
       fi 
    done
