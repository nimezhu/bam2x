for i in *.haplotype.*
    do
        j=${i%.*}
        k=${j%.*}
    mv $i $k
    done
