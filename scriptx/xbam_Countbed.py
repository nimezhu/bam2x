#!/usr/bin/python
# Wei Guifeng, <guifengwei@gmail.com>

import sys,os,argparse,time
import pysam


class Bed:
    def __init__(self,x):
        if isinstance(x,basestring):
            x=x.rstrip("\n\r").split("\t")
        try:
            x[-1].rstrip("\n\r")
        except:
            pass
        self.chr=x[0].strip()
        self.start=int(x[1])
        if self.start<0:
            self.start=0
        self.stop=int(x[2])
        try:
            self.id=x[3]
        except:
            self.id="NONAME"
        try:
            self.score=float(x[4])
        except:
            self.score=1
        try:
            self.strand=x[5]
        except:
            self.strand="."
        try:
            self.description="\t".join(x[6:])
        except:
            self.description=None
    def __str__(self):
        string=self.chr+"\t"+str(self.start)+"\t"+str(self.stop)+"\t"+str(self.id)+"\t"+str(self.score)
#       if(self.strand != "."):
#           string+="\t"+self.strand
        string+="\t"+str(self.strand)
        return string
    def length(self):
        return self.stop-self.start
    def overlap(A,B):
        if A is None or B is None : return 0
        if(A.chr != B.chr) : return 0
        if (A.stop <= B.start) : return 0
        if (B.stop <= A.start) : return 0
        return 1
    overlap=staticmethod(overlap)
    def merge(A,B):
        if A is None and B is None: return None
        if A and B: return Bed([A.chr,min(A.start,B.start),max(A.stop,B.stop),A.id+","+B.id,0,A.strand])
        if A: return A
        if B: return B
    merge=staticmethod(merge)
    def distance(self,bed):
        if(self.chr != bed.chr): return 10000000000
        if(Bed.overlap(self,bed)): return 0
        a=[abs(self.start-bed.stop),abs(self.stop-bed.stop),abs(self.start-bed.start),abs(self.stop-bed.start)]
        i=a[0]
        for j in a[1:]:
            if i>j:
                i=j
        return i
    def __cmp__(self,other):
        return cmp(self.chr,other.chr) or cmp(self.start,other.start) or cmp(self.stop,other.stop) or cmp(self.strand,other.strand)
    def getSeq(self,fn="/home/wyf/ws190/ce6.2bit"):
        if(fn is None): 
            print >>sys.stderr,"2Bit file not specified."
            return ""
        if self.start-self.stop>=0: return ""
        s=getSeq(fn,self.chr,self.start,self.stop)
        s=s.upper()
        if "-" in self.strand:
            s=Utils.rc(s)
        return s
    def upstream(self,bp):
        '''return the $bp bp upstream Bed Class Object'''
        chr=self.chr
        strand=self.strand
        id=self.id+"_up"+str(bp)
        if(self.strand=="+"):
            start=self.start-bp
            stop=self.start
        else:
            start=self.stop
            stop=self.stop+bp
        if (start<0):start=0
        x=[chr,start,stop,id,0,strand]
        return Bed(x)
    def downstream(self,bp):
        '''return the $bp bp downstream Bed Class Object'''
        chr=self.chr
        strand=self.strand
        id=self.id+"_down"+str(bp)
        if(self.strand=="+"):
            start=self.stop
            stop=self.stop+bp
        else:
            start=self.start-bp
            stop=self.start
        x=[chr,start,stop,id,0,strand]
        return Bed(x)
    def upstreamextend(self,bp):
        '''return the bed itself and $bp bp upstream Bed Class Object'''
        chr=self.chr
        strand=self.strand
        id=self.id+"_upextend"+str(bp)
        if(self.strand!="-"):
            start=self.start-bp
            stop=self.stop
        else:
            start=self.start
            stop=self.stop+bp
        x=[chr,start,stop,id,0,strand]
        return Bed(x)
    def downstreamextend(self,bp):
        '''return the bed itself and $bp bp downstream Bed Class Object'''
        chr=self.chr
        strand=self.strand
        id=self.id+"_downextend"+str(bp)
        if(self.strand!="-"):
            start=self.start
            stop=self.stop+bp
        else:
            start=self.start-bp
            stop=self.stop
        x=[chr,start,stop,id,0,strand]
        return Bed(x)
    def updownextend(self,bp):
        '''return the bed itself and $bp bp flanking regions Bed Class Object
'''
        chr=self.chr
        strand=self.strand
        id=self.id+"_updownextend"+str(bp)
        start=self.start-bp
        stop=self.stop+bp
        x=[chr,start,stop,id,0,strand]
        return Bed(x)
    def string_graph(self,scale=1):
        n=self.length()*scale
        n=int(n)
        s=""
        c="|"
        if(self.strand == "+"):
            c=">"
        if(self.strand == "-"):
            c="<"
        if (n==0): n=1
        for i in range(n): s+=c
        return s
    def tss(self):
        '''return the TSS class of transcription start site, name is geneid_tss'''
        pos=self.stop
        if self.strand=="+":
            pos=self.start
        else:
            pos=self.stop  #should minus one or not?
        return TSS([self.chr,pos,self.strand,self.id+"_tss"])
    def tts(self):
        if self.strand=="+":
            pos= self.stop
        else:
            pos=self.start
        return TSS([self.chr,pos,self.strand,self.id+"_tts"])
    def strand_cmp(self,bed):
        if bed.strand == ".":
            return "."
        if self.strand == bed.strand:
            return "+"
        else:
            return "-"
    def distance_to_tss(self,bed):
        if bed.chr != self.chr :
            return None
        tss=self.tss().tss
        pos=(bed.start+bed.stop)/2
        # strand=self.strand_cmp(bed)
        if(self.strand=="+"):
            rpos=pos-tss
        else:
            rpos=tss-pos
        return rpos
    def distance_to_tts(self,bed):
        if bed.chr != self.chr:
            return None
        tts=self.tts().tss
        pos=(bed.start+bed.stop)/2
        if(self.strand=="+"):
            rpos=pos-tts
        else:
            rpos=tts-pos
        return rpos




def parse_argument():
    ''' argument parser'''
    p=argparse.ArgumentParser(description='example: %(prog)s --bed *.bed --bam bamfile > output.tab',
                            epilog="dependency pysam")
    p.add_argument('--bed',dest='bed',metavar='Bed',type=str,
                    help="BED file input. Your interested genomic region")
    p.add_argument('--bam',dest='bam',type=str,nargs='+',
                    help=" bamfile or bamlist. eg. bamfile or 1.bam 2.bam ...")
    if len(sys.argv) == 1 :
        sys.exit(p.print_help())
    args=p.parse_args()
    return args

def CountReads(fn):
    """ CountReads in bamfile and return the total number of the aligned reads"""
    print >>sys.stderr,"counting ",fn
    samfile=pysam.Samfile(fn,"rb")
    ss=0
    for chr in samfile.references:
        print >>sys.stderr,"counting ",chr,"\r",
        ss += samfile.count(chr)
    samfile.close()
    return ss


def get_bamfile_readsnum(bamfile):
    '''INPUT: args.bamfile or args.RNA
       OUTPUT: bamlist and File_reads_num'''
    bamlist=[]
    File_reads_num={}
    if bamfile:
        for each in bamfile:
            if os.path.getsize(str(each)) < 1000000:
                f=open(each,'r')
                for line in f:
                    if line.startswith('#'): continue
                    line=line.strip()
                    a = line.split('\t')
                    bamlist.append(a[0])
                    if len(a)>1:
                        File_reads_num[a[0]]=long(a[1])
                    else:
                        File_reads_num[a[0]]=CountReads(a[0])
            else:
                f=open(each,'rb')
                bamlist.append(each)
                File_reads_num[each]=CountReads(each)
    return bamlist,File_reads_num

def read_bed(bedfile):
    ''' read the bed file'''
    bedlist=[]
    f=open(bedfile,'r')
    for line in f:
        if not line.startswith("#"):
            line=line.strip().split()
            a= Bed(line)
            bedlist.append(a)
    return bedlist

def BamInBed(bamfilename):
    '''
        Reading bamfile into tables, BED region
    '''
    print >>sys.stderr,"\nProcessing BAM File:",bamfilename,"into bed"
    samfile=pysam.Samfile(bamfilename,"rb")
    reads_number=File_reads_num[bamfilename]
    j=0
    for i in bedlist:
        j+=1
        if j%1000==0: print >>sys.stderr,j,"/",len(bedlist),"bed regions\r",
        c=0
        id=i.id
        if not expression.has_key(id):
            expression[id]=[]
        try:
            c=samfile.count(i.chr,i.start,i.stop)
        except:
            c=0
        rpkm=1000000000.0 *c/float(i.stop-i.start)/reads_number
        expression[id].append(rpkm)
    samfile.close()
    print >>sys.stderr," Done!"

def output_table():
    print "# Bedfile: ", args.bed
    for i,x in enumerate(bamlist):
        print "# BAM."+str(i+1)+":",x,"\t",File_reads_num[x],"reads"
    print "# chr\tstart\tstop\tid\tscore\tstrand",
    for i in range(len(bamlist)):
        print "\tBAM."+str(i+1),
    print
    for i,bed in enumerate(bedlist):
        print bed,
        for j in expression[bed.id]:
            print "\t%.3f"%j,
        print
    

def main():
    global args,File_reads_num,bedlist,bamlist,expression
    time
    File_reads_num={}
    expression={}

    ISOTIMEFORMAT='%Y-%m-%d %X'
    print "# xbam_countbed.py 0.1 : ",time.strftime( ISOTIMEFORMAT, time.localtime() )

    # parse the argument
    args=parse_argument()
    bamfile = args.bam
    bedfile = args.bed

    bedlist=read_bed(bedfile)
    bamlist,File_reads_num = get_bamfile_readsnum(bamfile)

    for i in bamlist:
        BamInBed(i)

    # output
    output_table()


if __name__=="__main__":
    main()


