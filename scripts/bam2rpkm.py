#!/usr/bin/env python
# programmer : zhuxp
# usage:
import sys
from getopt import getopt
import pysam


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
            pos=self.stop-1  #should minus one or not?
        #return TSS([self.chr,pos,self.strand,self.id+"_tss"])
        return TSS([self.chr,pos,self.strand,self.id])
    def tts(self):
        if self.strand=="+":
            pos= self.stop-1
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
    print >>sys.stderr,"bam2rpkm.py:  count reads in bed region"
    print >>sys.stderr,"Version:"+Version+"\n"
    print >>sys.stderr,"Library dependency: pysam\n\n"
    print >>sys.stderr,"Usage: bam2rpkm.py <options> --bed file.bed file1.bam file2.bam > output.tab"
    print >>sys.stderr,"       bam2rpkm.py <options> --bed file.bed --bam bamfiles > output.tab"
    print >>sys.stderr,"Example:  bam2rpkm.py --bed file.bed  file.bam > output.tab"
    print >>sys.stderr,"          bam2rpkm.py --bam bamlistfile --bed ncRNA.bed > output.tab"
    print >>sys.stderr,"                       bamlistfile example 1:"
    print >>sys.stderr,"                             /data/user/1.bam  10000000  # filename<tab>reads number<enter>"
    print >>sys.stderr,"                             /data/user/2.bam  11001100"
    print >>sys.stderr,"                       bamlistfile example 2:"
    print >>sys.stderr,"                             /data/user/1.bam"
    print >>sys.stderr,"                             /data/user/2.bam"
    print >>sys.stderr,"Options:"
    print >>sys.stderr,"   -h,--help          show this help message"
    print >>sys.stderr,"   -b binsize         TSS up and down bp number default: 2000  (-1000,+1000)"
    print >>sys.stderr,"   --bam bamlistfile"
    print >>sys.stderr,"                      a file contains bam filenames and reads number<optional> in each line"
    print >>sys.stderr,"   -r,--reads         not normalize to RPKM"
    print >>sys.stderr,"   "

def BamToRPKM(bamfilename,norm=1000000):
    print >>sys.stderr,"\nProcessing File:",bamfilename
    samfile=pysam.Samfile(bamfilename,"rb")
    reads_number=File_reads_number[bamfilename]
    j=0
    for i in bedList:
        j+=1
        if j%1000==0: print >>sys.stderr,j,"/",len(bedList),"bed entries","\r",
        try:
            c=samfile.count(i.chr,i.start,i.stop)
        except:
            c=0
        if ifNormalize:
            c=float(c)*norm*1000/reads_number/(i.stop-i.start)
        if not Table.has_key(i.id):
            Table[str(i.id)]=[]
        Table[i.id].append(c) 
    samfile.close()

def OutputTable():
    
    print "# Bamlist:"
    for i,x in enumerate(bamlist):
        print "#    BAM."+str(i+1)+":",x,"\t",File_reads_number[x],"reads"
    print "# Normalized to RPKM:",ifNormalize
    print "# chr\tstart\tend\tid\tscore\tstrand",
    for i in range(len(bamlist)):
        print "\tBAM."+str(i+1),
    print 
    print "# chr\tstart\tend\tid\tscore\tstrand",
    for i in range(len(bamlist)):
        dirfile=bamlist[i].split("/")
        file=dirfile[-1].replace(".bam","")
        print "\t",file,
    print 
    for i in bedList:
        print i,
        for j in Table[i.id]:
            if not ifNormalize:
                print "\t",j,
            else:
                print "\t%.3f"%j,
        print




def Main():
    global Table,Version,ifNormalize,File_reads_number,bedList,bamlist,Norm
    Version="0.1"
    Table={}
    if len(sys.argv)<2:
        show_help()
        exit(0)
    opts,restlist = getopt(sys.argv[1:],"hb:g:ru",\
                        ["help","bam=","reads","bed="])
    ifNormalize=True
    File_reads_number={}
    ifBamFileList=0
    bedList=[]
    Norm=1000000
    ifUniq=False
    BEDINDEX=0


    for o, a in opts:                       
        if o in ("-h","--help"): 
            show_help()
            exit
        if o in ("-b","--binsize"):
                Binsize=int(a)
        if o in ("--bam"):
            ifBamFileList=1
            BamFileList=a
        if o in ("--reads","-r"):
            ifNormalize=False
        if o in ("--bed"):
            bedfile=a
        if o in ("--uniq","-u"):
            ifUniq=True
    bedNames={};
    if 'bedfile' in dir():
        try:
            f=open(bedfile)
            print "# Bed File:",bedfile
            for line in f:
                line=line.strip()
                if line[0]=='#': continue
                x=line.split("\t")
                x=Bed(x)
                if bedNames.has_key(x.id):
                    BEDINDEX+=1
                    x.id=x.id+".IDX"+str(BEDINDEX)
                    bedList.append(x)
                else:
                    bedList.append(x)
                    bedNames[x.id]=1
                    BEDINDEX+=1
            f.close()
        except:
            print >>sys.stderr,"can't open bedfile",bedfile
            exit(0)
    bedList.sort()
    if len(bedList)==0:
        print >>sys.stderr,"No bed file assigned."
        exit(0)
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
                    print >>sys.stderr,"counting ",chr
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
    for i in bamlist:
        BamToRPKM(i,Norm)
    OutputTable()

    


    
if __name__=="__main__":
    Main()

