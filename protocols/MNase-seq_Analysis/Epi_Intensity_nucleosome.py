import sys,argparse,os
from xplib import TableIO
from xplib import DBI
from xplib.Annotation import Bed
import numpy as np

def ParseArg():
    p=argparse.ArgumentParser(description="get intensity of epi-modification on single nucleosome",epilog="library dependency: xplib")
    p.add_argument("-N","--Nucleosome",dest="nucleosome",type=str,help="xls file containing nucleosome location information from Danpos output")
    p.add_argument("-b","--beds",nargs='+',dest="beds",type=str,help="bed files for epigenetic data")
    p.add_argument("-l","--length",dest="len",type=int,default=200,help="average length of ChIP-seq fragment,default:200")
    p.add_argument("-n","--name",nargs='+',dest='name',type=str,help='name of each bed sample (to be wrote on the header)')
    p.add_argument("-o","--output",dest="output",type=str,help="output file name (can be .txt)")
    if len(sys.argv)==1:
        print >>sys.stderr,p.print_help()
        exit(0)
    return p.parse_args()

def find_center(read,half_len):

    # read.id in strand for bed with for columns
    if read.id=="+":
        center=read.start+half_len
    elif read.id=="-":
        center=read.stop-half_len
    else:
        print >>sys.stderr,read.id,"No strand"
    return center


'''
      mi
    -----
    |    \
    |     \
    |      \
    |       \
    --------
  Mid  ma

'''

global mi,ma
mi=75
ma=125


def Main():
    args=ParseArg()

    #store bed files with indexing and count information:
    bed={}

    print >>sys.stderr,"Starting index bed files:"
    for i in range(len(args.beds)):
        temp_name=args.name[i]
        print >>sys.stderr,"  #Indexing for bed file of",temp_name,"\r",
        bed[temp_name]=DBI.init(args.beds[i],'bed')

    half_len=int(args.len)
    print >>sys.stderr
    print >>sys.stderr,"Reading nucleosome peak xls file from Danpos."
    nucleosomes=TableIO.parse(args.nucleosome,'metabed',header=True)

    print >>sys.stderr,"Start Counting..."
    count_matrix=[]


    out=open(args.output,"w")
    line_head=open(args.nucleosome,'r').readline().strip()
    line_head=line_head+"\t"+"\t".join(str(f) for f in args.name)
    print >>out,line_head
    for i in nucleosomes:
        chrom=i.chr
      
        if chrom == 'chrY' or chrom == 'chrX' or chrom == 'chrM':
            continue
        center=int(i.start+i.end)/2
        count=np.zeros(len(args.beds),dtype="float")
        line=str(i)
        for k,name in enumerate(bed.keys()):
            for j in bed[name].query(Bed([chrom,center-ma-(half_len-75),center+ma+(half_len-75)])):
                j_center=find_center(j,half_len)
                weight = max(min(1,(ma-abs(j_center-center))/25.0),0)
                count[k]+=weight
        line = line + "\t" + "\t".join(str(f) for f in count)
        print >>out,line
        count_matrix.append(count)

if __name__=="__main__":
    Main() 


        
                        

'''
TableIO.parse("Pn_E14_merge_mm9.sorted.Fnor.ajClonal.smooth.peaks.xls",'metabed',header=True)
'''



