MNase-seq_Analysis:
-------------------

###  1. distogram and phasogram. ###
     Based on the idea of supplemental figure 3 of this paper:  
     http://www.nature.com/nature/journal/v474/n7352/full/nature10002.html  
     
     * Example command:
       * see run-phasogram-E14.sh
       * similar for distogram
     * Example output:
       * see pdf files


###  2. midpoint-gram ###
     Originated from Figure 1C of this paper:  
     http://www.plosgenetics.org/article/info%3Adoi%2F10.1371%2Fjournal.pgen.1003036  

     * Example command:
       ```
       python distogram_phasogram.py MNase_E14_mm9_rmdup.sorted.bam -n 50000000 -m -x -200 200 -i 1 -o E14_rmdup.pdf
       ```
     * Example output:
       * see pdf files


###  3. Count epi- mark data in single nucleosome level ###
     First achieve nucleosome position information using Danpos software.
     Danpos: https://code.google.com/p/danpos/
     The script take output xls file of Danpos as input file for nucleosome positions.
     The figure file (histone mark intensity-single nucleosome.jpg) illustrate the idea.

     Example command:
       see Epi_Intensity_nucleosome_example.sh



## Program argument options: ##

  #### 1. distogram_phasogram.py ####
```
usage: distogram_phasogram.py [-h] [-d | -p | -m] [-b input_bam] [-o OUTPUT]
                              [-i PILE] [-x XLIM [XLIM ...]] [-n NUM]
                              input_rmdup

draw distogram/phasogram/Midpoint for MNase_seq data

positional arguments:
  input_rmdup           dup removed input bam file for mapped MNase-seq data

optional arguments:
  -h, --help            show this help message and exit
  -d, --distogram       draw distogram for MNase_seq data
  -p, --phasogram       draw phasogram for MNase-seq data
  -m, --midpoint        draw the distribution of MNase-seq midpoint
  -b input_bam, --input_o input_bam
                        original input bam file for mapped MNase-seq
                        data,necessary when pile>1
  -o OUTPUT, --output OUTPUT
                        the output figure file, can be format of emf, eps,
                        pdf, png, ps, raw, rgba, svg, svgz
  -i PILE, -pile PILE   conditioning the analysis on sites with #[pile] or
                        more read starts (default=1)
  -x XLIM [XLIM ...], --xlim XLIM [XLIM ...]
                        range for x axis: min_x, the left bound; max_x, the
                        right bound. (default: -200,2000)
  -n NUM, --num NUM     number of distances included to draw density function
                        (default: 5000000)

Library dependency: pysam, matplotlib, scipy, numpy
```
  
####  2. plot_histogram.py ####
```
usage: plot_histogram.py [-h] [-d | -p] [-x XLIM [XLIM ...]] [-o OUTPUT]
                         input_count

draw disto/phasogram for processed count data

positional arguments:
  input_count           input txt file for processed count data

optional arguments:
  -h, --help            show this help message and exit
  -d, --distogram       draw distogram for MNase_seq data
  -p, --phasogram       draw phasogram for MNase-seq data
  -x XLIM [XLIM ...], --xlim XLIM [XLIM ...]
                        range for x axis: min_x, the left bound; max_x, the
                        right bound. (default: -200,2000)
  -o OUTPUT, --output OUTPUT
                        the output figure file format, can be format of emf,
                        eps, pdf, png, ps, raw, rgba, svg, svgz

Library dependency: matplotlib, scipy, numpy


  3. Epi_Intensity_nucleosome.py
usage: Epi_Intensity_nucleosome.py [-h] [-N NUCLEOSOME] [-b BEDS [BEDS ...]]
                                   [-l LEN] [-n NAME [NAME ...]] [-o OUTPUT]

get intensity of epi-modification on single nucleosome

optional arguments:
  -h, --help            show this help message and exit
  -N NUCLEOSOME, --Nucleosome NUCLEOSOME
                        xls file containing nucleosome location information
                        from Danpos output
  -b BEDS [BEDS ...], --beds BEDS [BEDS ...]
                        bed files for epigenetic data
  -l LEN, --length LEN  average length of ChIP-seq fragment,default:300
  -n NAME [NAME ...], --name NAME [NAME ...]
                        name of each bed sample (to be wrote on the header)
  -o OUTPUT, --output OUTPUT
                        output file name (can be .txt)

library dependency: xplib
```
