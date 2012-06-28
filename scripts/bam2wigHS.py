#!/usr/bin/python
# Programmer : zhuxp
# Date: 
# Last-modified: 28 Jun 2012 14:51:25
import os,sys,argparse
import pysam
def ParseArg():
    ''' This Function Parse the Argument '''
    p=argparse.ArgumentParser( description = 'Example: %(prog)s -i file.bam -o file.wig', epilog='Library dependency : pysam')
    p.add_argument('-v','--version',action='version',version='%(prog)s 0.1')
 
 
    p.add_argument('-n','--normalize',dest="Normalize",action="store_true",help="normalize to 100 bp * 10_000_000 reads")
    p.add_argument('--input','-i',type=str,dest='Bamfile',help="the input alignment bamfile")
    p.add_argument('--output','-o',type=str,dest='Output',default="stdout", help="the output wig file  default: %(default)s")
    if len(sys.argv)==1:
        print >>sys.stderr,p.print_help()
        exit(0)

    return p.parse_args()
def ArrayToWig(a):
    state=0
    start=0
    stop=0
    thres=100
    zero_count=0
    segment=[]
    for i,x in enumerate(a):
        if x>0 and state==0: 
            start=i
            state=1
            zero_count=0
        if x>0 and state==1:
            stop=i
            zero_count=0
        if x==0 and state==0:
            continue
        if x==0 and state==1:
            zero_count+=1
            if zero_count > thres:
                state=0
                segment.append((start,stop))
    if state==1:
        segment.append((start,stop))
    return segment









    

def Main():
    global args,Table
    args=ParseArg()
    if args.Output=="stdout":
        out=sys.stdout
    else:
        out=open(args.Output,"w")
    samfile=pysam.Samfile(args.Bamfile,"rb")
    chrs=samfile.references
    lengths=samfile.lengths
    Table=[]
    print >>sys.stderr,"initialize table in memory"
    for i,x in enumerate(chrs):
        l=lengths[i]
        bins=l>>20
        Table.append([])  #Table[chr_id]
        for j in range(bins+1):
            Table[i].append([])
    print >>sys.stderr,"initial table done                                  "
    unmap=0
    map=0
    total_map_bps=0
    t=0
    print >>sys.stderr,"reading file",args.Bamfile,"into table"

    for sam in samfile:
        t+=1
        if t%100000==0: print >>sys.stderr,t,"reads\r",
        if sam.tid<0:
            unmap+=1
            continue
        i=sam.tid
        s=sam.pos>>20
        if len(Table[i][s])==0:
            Table[i][s]=[0 for row in range(1<<20)]
        e=(sam.pos+sam.qlen)>>20
        if e>s:
            Table[i][e]=[0 for row in range(1<<20)]
        for j in range(sam.pos,sam.pos+sam.qlen):
            k=j
            k=k>>20
            l=j-(k<<20)
            Table[i][k][l]+=1
        total_map_bps+=sam.qlen
        map+=1
    print >>sys.stderr,"Printing Wig\r",
    print >>out,"track type=wiggle_0 name="+args.Bamfile+" description="+args.Bamfile
    for i,chr in enumerate(chrs):
        for j in range(len(Table[i])):
            if len(Table[i][j])>0:
                S=ArrayToWig(Table[i][j])
                for a in S:
                    (start,stop)=a
                    print >>out,"fixedStep\tchrom="+chr,"start="+str(start),"step=1"
                    for k in range(start,stop+1):
                        print >>out,Table[i][j][k]


    print >>sys.stderr,"Printing Wig Done"


    
if __name__=="__main__":
    Main()


