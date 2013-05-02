#/usr/bin/env sh

if [ $# -ge 1 ]
    then
      dir=$1
    else
      dir="./"
fi
for i0 in $dir/*.fq
    do
       # dir=${i0%/*.fq}
        
        i=${i0##*/}
        if [ "$i" == "*.fq" ]
            then
              echo "no fq in this directory"
              break
        fi
        j2=${i%_1.fq}
        mate2="${j2}_2.fq"
        j1=${i%_2.fq}
        mate1="${j1}_1.fq"
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
