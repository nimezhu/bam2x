#!/usr/bin/python
# programmer : zhuxp
# usage:
import sys
from getopt import getopt
import pysam
def show_help():
    print >>sys.stderr,"bam2wig.py:  convert bamfile into wigfile"
    print >>sys.stderr,"Version:"+Version+"\n"
    print >>sys.stderr,"Usage: bam2wig.py <options> file.bam > file.wig"
    print >>sys.stderr,"Example:  bam2wig.py -b binsize file.bam > file.wig"
    print >>sys.stderr,"          bam2wig.py --bam bamlistfile --region chr1:1-10000 > file.wig"
    print >>sys.stderr,"                       bamlistfile example 1:"
    print >>sys.stderr,"                             /data/user/1.bam  10000000  # filename<tab>reads number<enter>"
    print >>sys.stderr,"                             /data/user/2.bam  11001100"
    print >>sys.stderr,"                       bamlistfile example 2:"
    print >>sys.stderr,"                             /data/user/1.bam"
    print >>sys.stderr,"                             /data/user/2.bam"
    print >>sys.stderr,"          bam2wig.py --bam bamlistfile --output_bamlist new_bamlistfile"
    print >>sys.stderr,"                       count reads number in bams"
    print >>sys.stderr,"Library dependency: pysam\n\n"
    print >>sys.stderr,"Options:"
    print >>sys.stderr,"   -h,--help          show this help message"
    print >>sys.stderr,"   -b binsize         window size for each bin, default: 200"
    print >>sys.stderr,"   --bam bamlistfile"
    print >>sys.stderr,"                      a file contains bam filenames and reads number<optional> in each line"
    print >>sys.stderr,"   --output_bamlist bamlist"
    print >>sys.stderr,"                      output a file contains bam filenames and reads number"
    print >>sys.stderr,"                      if this option is selected , do not output wig."
    print >>sys.stderr,"   --region chr:start-stop"
    print >>sys.stderr,"                      generate wig in selected region"
    print >>sys.stderr,"   --normalize        normalize to RPKM"
    print >>sys.stderr,"   "

def BamToRegionWig(bamfilename,region,norm=1000000):
    a=region.split(":")
    chr=a[0].strip()
    b=a[1].split("-")
    start=long(b[0])
    stop=long(b[1])
    
    samfile=pysam.Samfile(bamfilename,"rb")
    print >>sys.stderr,"Processing ",bamfilename
    reads_number=File_reads_number[bamfilename]
    lengths=samfile.lengths
    s=0
    for i in lengths:
        s+=i
    total_bins_num=s/Binsize
    l=float(reads_number)/total_bins_num

    bins_number=(stop-start)/Binsize
    if (stop-start)/Binsize!=0: bins_number+=1
    bins=[0 for row in range(bins_number)]

    sams=samfile.fetch(chr,start,stop)
    for sam in sams:
        pos=sam.pos+sam.qlen/2
        bin=(pos-start)/Binsize
        if bin>0 and bin<len(bins):
            bins[bin]+=1
    print "# Bamfile: ",bamfilename
    print "# Total Number Reads:",reads_number
    if not ifNormalize:
        print "# Lambda: %.4f"%l
    else:
        print "# Normalize To RPKM: True"
    a=bamfilename.split("/")
    bamfn=a[-1]
    bamfn=bamfn.replace("wgEncode","")
    bamfn=bamfn.replace("Histone","")
    bamfn=bamfn.replace("Broad","")
    bamfn=bamfn.replace("Sydh","")
    bamfn=bamfn.replace(".bam","")
    print "track type=wiggle_0 name=\""+bamfn+"\" description=\""+bamfn+"\" visibility=full autoScale=on color=150,150,255" 
    print "variableStep chrom="+chr+" span="+str(Binsize)
    for i in range(len(bins)):
        if bins[i]!=0:
            if ifNormalize:
                xx=float(bins[i])*norm*1000/reads_number/Binsize
                print start+i*Binsize,"\t%.4f"%xx
            else:
                print start+i*Binsize,"\t",bins[i]


def BamToWig(bamfilename,norm=1000000):
    hChr={}
    samfile=pysam.Samfile(bamfilename,"rb")
    a=bamfilename.split("/")
    bamfn=a[-1]
    print >>sys.stderr,"reading file",bamfilename,"now"
    chrs=samfile.references
    lengths=samfile.lengths
    total_bins_num=0
    for i in range(len(chrs)):
        bins=lengths[i]/Binsize
        if lengths[i]%Binsize!=0: bins+=1
        total_bins_num+=bins
        hChr[chrs[i]]=[0 for row in range(bins)]
    i==0
    for s in samfile:
        pos=s.pos+s.qlen/2
        hChr[chrs[s.tid]][pos/Binsize]+=1
        i+=1
        if i%100000==0: print >>sys.stderr,"\treading ",i,"reads now"
    print "# Bamfile: ",bamfilename
    print "# Total Reads Number:",i
    total_reads_num=i
    l=float(i)/total_bins_num
    if not ifNormalize: 
        print "# Lambda: %.4f"%l
    else:
        print "# Normalize to RPKM"
    print "track type=wiggle_0 name=\""+bamfn+"\" description=\""+bamfn+"\" visibility=full autoScale=on color=50,150,255" 
    for i in range(len(chrs)):
        print "variableStep chrom="+chrs[i]+" span="+str(Binsize)
        for j,x in enumerate(hChr[chrs[i]]):
            if x!=0:
                if ifNormalize:
                    xx=float(x)*norm*1000/total_reads_num/Binsize
                    print j*Binsize,"\t%.4f"%xx
                else:
                    print j*Binsize,"\t",x

def Main():
    global Version,Binsize,ifRegion,Region,ifNormalize,File_reads_number
    Version="0.1"
    if len(sys.argv)<2:
        show_help()
        exit
    opts,restlist = getopt(sys.argv[1:],"hb:r:n",\
                        ["help","binsize=","region=","bam=","normalize","output_bamlist="])
    Binsize=200
    ifRegion=0
    ifNormalize=0
    ifOutputBamlist=0
    File_reads_number={}
    ifBamFileList=0


    for o, a in opts:                       
        if o in ("-h","--help"): 
            show_help()
            exit
        if o in ("-b","--binsize"):
                Binsize=int(a)
        if o in ("-r","--region"):
                ifRegion=1
                Region=a
        if o in ("--bam"):
            ifBamFileList=1
            BamFileList=a
        if o in ("-n","--normalize"):
            ifNormalize=1
        if o in ("--output_bamlist"):
            output_bamlist_file=a
            ifOutputBamlist=1

    bamlist=[]
    if ifBamFileList:
        f=open(BamFileList)
        for i in f:
            i=i.strip()
            a=i.split("\t")
            bamlist.append(a[0])
            if len(a)>1:
                File_reads_number[a[0]]=long(a[1])   # using total reads number pre-computed in bamlistfile
            else:
                print >>sys.stderr,"counting ",a[0]  # counting reads number in bam , might be slow. 
                samfile=pysam.Samfile(a[0],"rb")
                ss=0
                for chr in samfile.references:
                    ss+=samfile.count(chr)
                File_reads_number[a[0]]=ss
                samfile.close()
    for i in restlist:
        i=i.strip()
        a=i.split('.')
        if a[-1]=="bam": 
            bamlist.append(i)
        samfile=pysam.Samfile(i,"rb")
        ss=0
        print >>sys.stderr,"counting ",i
        for chr in samfile.references:
            ss+=samfile.count(chr)
        File_reads_number[i]=ss
        samfile.close()
    if ifOutputBamlist:
        f=open(output_bamlist_file,"w")
        for i in bamlist:
            f.write(i+"\t"+str(File_reads_number[i])+"\n")
    else: 
        if ifRegion:
            print "browser position",Region
            print "browser hide all"
        for i in bamlist:
            if ifRegion:
                BamToRegionWig(i,Region)
            else:
                BamToWig(i)
    


    
if __name__=="__main__":
    Main()

