#!/usr/bin/python
# programmer : zhuxp
# usage:
import sys
from getopt import getopt
import pysam
from scipy.stats.stats import ttest_ind



class Bed:
    def __init__(self,x):
        
	self.chr=x[0].strip()
	self.start=int(x[1])
	self.stop=int(x[2])
        try:
	    self.id=x[3].strip()
        except:
            self.id="NONAME"
        try:
            self.score=float(x[4])
        except:
            self.score=0
	try:
	    self.strand=x[5].strip()
	except:
	    self.strand="."
    def __str__(self):
	string=self.chr+"\t"+str(self.start)+"\t"+str(self.stop)+"\t"+str(self.id)+"\t"+str(self.score)
#	if(self.strand != "."):
#	    string+="\t"+self.strand
	string+="\t"+str(self.strand)
	return string
    def length(self):
        return self.stop-self.start
    def overlap(A,B):
	if(A.chr != B.chr) : return 0
	if (A.stop <= B.start) : return 0
	if (B.stop <= A.start) : return 0
	return 1
    overlap=staticmethod(overlap)
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
    def core_promoter(self,bp=1000,down=500):
        '''return the $bp bp upstream Bed Class Object'''
        chr=self.chr
        strand=self.strand
        id=self.id+"_core_promoter"
        if(self.strand=="+"):
            start=self.start-bp
            stop=self.start+down
        else:
            start=self.stop-down
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



def show_help():
    print >>sys.stderr,"bam2tab.py\tVersion:"+Version
    print >>sys.stderr,"Usage: bam2tab.py <options> file1.bam file2.bam ... annotation_file1.bed annotation_file2.bed ... > file.tab"
    print >>sys.stderr,"Example: bam2tab.py 1.bam 2.bam > file.tab"
    print >>sys.stderr,"         bam2tab.py --bam bamlistfile > file.tab\n\n"
    print >>sys.stderr,"Options:"
    print >>sys.stderr,"\t--bin binsize[default 200] "
    print >>sys.stderr,"\t\t\twindow size for each bin, default is 200 "
    print >>sys.stderr,"\t--bam bamlistfile"
    print >>sys.stderr,"\t\t\ta file contains bam filenames [load many bamfiles]"
    print >>sys.stderr,"\t--readsnumber"
    print >>sys.stderr,"\t\t\tnot normalized to RPKM, report the reads number instead [default: False]"
    print >>sys.stderr,"\t-n,--normalize"
    print >>sys.stderr,"\t\t\tnormalize to RPKM"
    print >>sys.stderr,"\t-s,--slim"
    print >>sys.stderr,"\t\t\tnot print all zero line, output will be smaller using this option"
    print >>sys.stderr,"\t-c N,--case N"
    print >>sys.stderr,"\t\t\tthe first N bams are cases, compare them with the rest of bamfiles as control using standard t-test"
    exit(0)

def count_bins():
    if ifRegion:  #
        pass
    else:
        a=chromeSizes
        s=0
        for i,x in enumerate(a[1]):
            Offset[a[0][i]]=s
            bins=x/Binsize
            if x%Binsize!=0: bins=bins+1
            s+=bins
        return s
def pos2bin(chrome,pos):
    return Offset[chrome]+pos/Binsize
def bin2pos(bin):
    a=chromeSizes
    for i,x in enumerate(a[0]):
        if bin < Offset[x]:
            i=i-1
            break
    local_bin=bin-Offset[a[0][i]]
    s=Bed([a[0][i],local_bin*Binsize,local_bin*Binsize+Binsize-1,bin,".","."])
    return s



    

def readIntoTable(i,bamfile):
    a=chromeSizes
    samfile=pysam.Samfile(bamfile,"rb")
    j=0
    no_map=0
    for sam in samfile:
        j+=1
        if j%1000000==0: 
            print >>sys.stderr," reading ",j,"reads","\r",
        tid=sam.tid
        if tid<0: 
            no_map+=1
            continue
        chrome=samfile.getrname(sam.tid)
        pos=sam.pos+len(sam.seq)/2
        bin_index=pos2bin(chrome,pos)
 #       print chrome,pos,bin_index
        Table[i][bin_index]+=1
        TotalReadsNumber[i]+=1
    print >>sys.stderr," Mapping Reads:", TotalReadsNumber[i]
    print >>sys.stderr," No Mapping Reads:", no_map
def readBedIntoTable(i,bedfile):
    try:
        f=open(bedfile)
    except:
        print >>sys.stderr,"can't open file ",bedfile
        exit(0)
    for line in f:
        x=Bed(line.split("\t"))
        if not Offset.has_key(x.chr): continue
        bin_indexes=bed2bins(x.chr,x.start,x.stop)
        for j in bin_indexes:
            AnnoTable[i][j]=1

def bed2bins(chr,start,stop):
    a=chromeSizes
    j=0
    for i,x in enumerate(a[0]):
        if x==chr: j=i
    if stop > a[1][j] : stop=a[1][j]
    return range(Offset[chr]+start/Binsize,Offset[chr]+stop/Binsize+1)

def normalizeTable(normalized_readsnum=1000000):
    for i in range(len(Table)):
        for j in range(len(Table[i])):

            Table[i][j]=float(Table[i][j])*normalized_readsnum*1000/TotalReadsNumber[i]/Binsize
    
def print_header():
    pass
def print_table():
    print >>sys.stderr,"printing table now"
    for i in range(len(Table[0])):
        a=bin2pos(i)
        if ifSlim:
            s=0
            for j in range(len(Table)):
                s=s+Table[j][i]
            if s==0:
                continue
        print a.chr,"\t",a.start,"\t",a.stop,
        if ifCompare:
            x=[]
            y=[]
            for j in range(case_number):
                x.append(Table[j][i])
            for j in range(case_number,len(Table)):
                y.append(Table[j][i])
            (t,p)=ttest_ind(x,y)
            print "\t",t,"\t",p,
        for j in range(len(Table)):
            if ifNormalize:
                print "\t%.2f"%Table[j][i], 
            else:
                print "\t",Table[j][i],
        for j in range(len(AnnoTable)):
            print "\t",AnnoTable[j][i],
        print 



def Main():
    global Version,Binsize,Table,ifRegion,Offset,TotalReadsNumber,AnnoTable,case_number,ifCompare,ifSlim,chromeSizes,ifNormalize
    Version="0.03"
    if len(sys.argv)==1: show_help()
    opts,restlist = getopt(sys.argv[1:],"ohr:nc:s",\
                        ["readsnumber","help","region=","bin=","normalize","gene=","case_number=","bam=","slim"])
    
    Offset={}
    Table=[]
    ifRegion=False
    ifNormalize=True
    Binsize=200
    normalized_readsnum=10000000
    AnnoTable=[]
    case_number=1
    ifCompare=0
    ifFbam=0
    ifSlim=0

    for o, a in opts:
        if o in ("-h","--help"): show_help()
       # if o in ("--gene"): Gene = a
       # if o in ("-r","--region"): 
       #     region=a
       #     ifRegion=True
        if o in ("--readsnumber"):
            ifNormalize=False
        if o in ("--bin"): 
            Binsize=int(a)
        if o in ("-n","--normalize"):
            ifNormalize=True
        if o in ("-c","--case_number"):
            case_number=int(a)
            ifCompare=1
        if o in ("--bam"):
            bamfiles_file=a
            ifFbam=1
        if o in ("--slim","-s"):
            ifSlim=1


    bam_filenames=[]
    bed_filenames=[]
    for i in restlist:
        a=i.split('.')
        if a[-1]=="bed": bed_filenames.append(i)
        if a[-1]=="bam": bam_filenames.append(i)
    if ifFbam:
        f=open(bamfiles_file)
        for line in f:
            line=line.strip()
            a=line.split("\t")
            bam_filenames.append(a[0])
    bam_num=len(bam_filenames)
    bed_num=len(bed_filenames)
    TotalReadsNumber=[0 for row in range(bam_num)]
#### parse chromosome size from bam file    
    for bamfile in bam_filenames:
        bamfile1=pysam.Samfile(bamfile,"rb")
        tempChromeSizes=(bamfile1.references,bamfile1.lengths)
        bamfile1.close()
        if not 'chromeSizes' in dir():
            chromeSizes=tempChromeSizes
        if len(chromeSizes[0])<len(tempChromeSizes[0]):
            chromeSizes=tempChromeSizes                   #in case that some bam don't have chromosome Y

#### parse end



    Table=[[0.0 for col in range(count_bins())] for row in range(bam_num)]
    AnnoTable=[[0 for col in range(count_bins())] for row in range(bed_num)]
    #  print_header()
    for i,x in enumerate(bed_filenames):
        print >>sys.stderr,"\nreading annotation file",x,"now"
        readBedIntoTable(i,x)
 
    for i,x in enumerate(bam_filenames):
        print >>sys.stderr,"\nreading file", x ,"now."
        readIntoTable(i,x)
    if ifNormalize:
        normalizeTable(normalized_readsnum)
    print "#  Binsize:",Binsize
    print "#  Normalized To RPKM:",ifNormalize
    print "#  BAM Files:",len(bam_filenames)
    for i,x in enumerate(bam_filenames):
      if ifCompare:
          if i<case_number:
              print "#\t\tCase",i+1,":",
          else:
              print "#\t\tControl",i+1-case_number,":",
      else:
          print "#\t\tNo,",i+1,":", 
      print x,"\t",TotalReadsNumber[i],"reads"
    print "#  BED Files:",len(bed_filenames)
    for i,x in enumerate(bed_filenames):
      print "#          No",i+1,":",x
    
    
    #################################
    print "#chr\tstart\tstop",
    if ifCompare:
        print "\tt-value\tp-value"
        for i in range(case_number):
            print "\tCase_"+str(i+1),
        for i in range(case_number,len(bam_filenames)):
            print "\tControl_"+str(i-case_number+1),

    else:
        for i in range(len(bam_filenames)):
            print "\tBam_"+str(i+1),
    for i in range(len(bed_filenames)):
        print "\tBed_"+str(i+1),
    print
    ##################################
  
    print_table()
    print >>sys.stderr,"Done!"







    
if __name__=="__main__":
    Main()

