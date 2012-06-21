#!/usr/bin/python
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
        self.uid=self.chr+":"+str(self.start)+"-"+str(self.stop)+"_"+self.id
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

        




class GeneBed(Bed):
    '''Gene Class'''
    def __init__(self,x):
        try:
            self.bin=int(x[0])
            if(self.bin < 100000):
                x=x[1:]
        except:
            pass
        self.id=x[0].rstrip()
        self.chr=x[1].rstrip()
        self.strand=x[2]
        if(self.strand == 1):
            self.strand="+"
        elif(self.strand == -1):
            self.strand="-"
        elif(self.strand == 0):
            self.strand="."
        self.start=int(x[3])
        self.stop=int(x[4])
        self.cds_start=int(x[5])
        self.cds_stop=int(x[6])
        self.exon_count=int(x[7])
        self.exon_starts=x[8].strip(",").split(",")
        self.exon_stops=x[9].strip(",").split(",")
        for i in range(self.exon_count):
            try:
                self.exon_starts[i]=int(self.exon_starts[i])
                self.exon_stops[i]=int(self.exon_stops[i])
            except:
                pass
        self.score=0
        try:
            self.protein_id=x[10]
            self.name2=x[10]
        except:
            self.name2="None"
        try:
            self.align_id=x[11]
        except:
            pass
        self.uid=self.chr+":"+str(self.start)+"-"+str(self.stop)+"_"+self.id
    def _exon(self,i):
        '''internal fucntion to call the exon position'''
        if i > self.exon_count:
            return None
        if self.strand=="+":
            return self.exon_starts[i-1],self.exon_stops[i-1]
        if self.strand=="-":
            n=self.exon_count
            return self.exon_starts[n-i],self.exon_stops[n-i]
    def getExon(self,i):
        '''get exon and return Bed Class of exon, the first exon is Exon_1'''
        if(i > self.exon_count or i < 1):return None
        start,stop=self._exon(i)
        id=self.id+"_"+"Exon_"+str(i)
        x=[self.chr,start,stop,id,0,self.strand]
        return Bed(x)
    def Exons(self):
        '''return a list of Bed classes of exon'''
        a=[]
        for i in range(self.exon_count):
            a.append(self.getExon(i+1))
        return a
    def cdna_length(self):
        '''sum of the length of exons'''
        s=0
        for i in self.Exons():
            s=s+i.length()
        return s
    def _intron(self,i):
        if i > self.exon_count-1:
            return None
        if self.strand=="+":
            return self.exon_stops[i-1],self.exon_starts[i]
        if self.strand=="-":
            n=self.exon_count
            return self.exon_stops[n-i-1],self.exon_starts[n-i]
    def getIntron(self,i):
        '''get exon and return Bed Class of intron, the first intron is Intron_1'''
        if(i>self.exon_count-1 or i< 1): return None
        start,stop=self._intron(i)
        id=self.id+"_"+"Intron_"+str(i)
        x=[self.chr,start,stop,id,0,self.strand]
        return Bed(x)
    def Introns(self):
        '''return a list of Bed classes of exon'''
        a=[]
        for i in range(self.exon_count-1):
            a.append(self.getIntron(i+1))
        return a
    def PrintString(self):
        print self.chr,
        print self.id,self.exon_starts,self.exon_stops
    def cds(self): 
        id=self.id+"_"+"cds"
        chr=self.chr
        start=self.cds_start
        stop=self.cds_stop
        strand=self.strand
        cds_start=self.cds_start
        cds_stop=self.cds_stop
        exon_count=1
        exon_starts=str(cds_start)+","
        exon_stops=""
        for i in range(self.exon_count):
            if (self.exon_starts[i] > self.cds_start and self.exon_starts[i] < self.cds_stop):
                    exon_starts+=str(self.exon_starts[i])+","
                    exon_count+=1
            if self.exon_stops[i] > self.cds_start and self.exon_stops[i] < self.cds_stop:
                    exon_stops+=str(self.exon_stops[i])+","
        exon_stops+=str(cds_stop)+","
        return GeneBed([id,chr,strand,start,stop,cds_start,cds_stop,exon_count,exon_starts,exon_stops])
    def utr5(self): 
        if(self.strand == "+"):
            if(self.cds_start==self.start):
                return None
            id=self.id+"_"+"utr5"
            chr=self.chr
            start=self.start
            stop=self.cds_start
            strand=self.strand
            cds_start=self.cds_start
            cds_stop=self.cds_start
            exon_count=0
            exon_stop_count=0
            exon_starts=""
            exon_stops=""
            for i in range(self.exon_count):
                if self.exon_starts[i] < self.cds_start:
                    exon_starts+=str(self.exon_starts[i])+","
                    exon_count+=1
                if self.exon_stops[i] < self.cds_start:
                    exon_stops+=str(self.exon_stops[i])+","
                    exon_stop_count+=1
            if exon_stop_count < exon_count:
                    exon_stops+=str(self.cds_start)+","
            return GeneBed([id,chr,strand,start,stop,cds_start,cds_stop,exon_count,exon_starts,exon_stops])
        if(self.strand == "-"):
            if(self.cds_stop==self.stop):
                return None
            id=self.id+"_"+"utr5"
            chr=self.chr
            start=self.cds_stop
            stop=self.stop
            strand=self.strand
            cds_start=self.cds_stop
            cds_stop=self.cds_stop
            exon_count=0
            exon_start_count=0
            exon_starts=""
            exon_stops=""
            for i in range(self.exon_count):
                if self.exon_starts[i] > self.cds_stop:
                    exon_starts+=str(self.exon_starts[i])+","
                    exon_start_count+=1
                if self.exon_stops[i] > self.cds_stop:
                    exon_stops+=str(self.exon_stops[i])+","
                    exon_count+=1
            if exon_start_count < exon_count:
                    exon_starts=str(self.cds_stop)+","+exon_starts
            return GeneBed([id,chr,strand,start,stop,cds_start,cds_stop,exon_count,exon_starts,exon_stops])
    def utr3(self): 
        if(self.strand == "-"):
            if(self.cds_start==self.start):
                return None
            id=self.id+"_"+"utr3"
            chr=self.chr
            start=self.start
            stop=self.cds_start
            strand=self.strand
            cds_start=self.cds_start
            cds_stop=self.cds_start
            exon_count=0
            exon_stop_count=0
            exon_starts=""
            exon_stops=""
            for i in range(self.exon_count):
                if self.exon_starts[i] < self.cds_start:
                    exon_starts+=str(self.exon_starts[i])+","
                    exon_count+=1
                if self.exon_stops[i] < self.cds_start:
                    exon_stops+=str(self.exon_stops[i])+","
                    exon_stop_count+=1
            if exon_stop_count < exon_count:
                    exon_stops+=str(self.cds_start)+","
            return GeneBed([id,chr,strand,start,stop,cds_start,cds_stop,exon_count,exon_starts,exon_stops])
        if(self.strand == "+"):
            if(self.cds_stop==self.stop):
                return None
            id=self.id+"_"+"utr3"
            chr=self.chr
            start=self.cds_stop
            stop=self.stop
            strand=self.strand
            cds_start=self.cds_stop
            cds_stop=self.cds_stop
            exon_count=0
            exon_start_count=0
            exon_starts=""
            exon_stops=""
            for i in range(self.exon_count):
                if self.exon_starts[i] > self.cds_stop:
                    exon_starts+=str(self.exon_starts[i])+","
                    exon_start_count+=1
                if self.exon_stops[i] > self.cds_stop:
                    exon_stops+=str(self.exon_stops[i])+","
                    exon_count+=1
            if exon_start_count < exon_count:
                    exon_starts=str(self.cds_stop)+","+exon_starts
            return GeneBed([id,chr,strand,start,stop,cds_start,cds_stop,exon_count,exon_starts,exon_stops])
    def toBedString(self):
        return Bed.__str__(self)
    def strand_cmp(self,bed):
        if bed.strand == ".":
            return "."
        if self.strand == bed.strand:
            return "+"
        else:
            return "-"
    def region_cmp(self,bed):
        '''report Bed Object overlap with which exon and intron or don't overlap with gene'''
        if not Bed.overlap(self,bed):
            r="Non-overlap"
            return r
        r="Overlap: "
        for i in range(self.exon_count):
            start,stop=self._exon(i+1)
            if(bed.start<stop and start<bed.stop):
                r+="Exon_"+str(i+1)+","
        for i in range(self.exon_count-1):
            start,stop=self._intron(i+1)
            if(bed.start < stop and start < bed.stop):
                r+="Intron_"+str(i+1)+","
        return r

                    
        
def show_help():
    print >>sys.stderr,"bamInBed.py:  count reads in bed region"
    print >>sys.stderr,"Version:"+Version+"\n"
    print >>sys.stderr,"Library dependency: pysam\n\n"
    print >>sys.stderr,"Usage: bamInBed.py <options> --gene gene.tab file1.bam file2.bam > output.tab"
    print >>sys.stderr,"       bamInBed.py <options> --bed gene.bed --bam bamfiles > output.tab"
    print >>sys.stderr,"Example:  bamInBed.py --gene hg19.refseq.tab -b binsize file.bam > output.tab"
    print >>sys.stderr,"          bamInBed.py --bam bamlistfile --bed ncRNA.bed > output.tab"
    print >>sys.stderr,"                       bamlistfile example 1:"
    print >>sys.stderr,"                             /data/user/1.bam  10000000  # filename<tab>reads number<enter>"
    print >>sys.stderr,"                             /data/user/2.bam  11001100"
    print >>sys.stderr,"                       bamlistfile example 2:"
    print >>sys.stderr,"                             /data/user/1.bam"
    print >>sys.stderr,"                             /data/user/2.bam"
    print >>sys.stderr,"Options:"
    print >>sys.stderr,"   -h,--help          show this help message"
    print >>sys.stderr,"   --bam bamlistfile"
    print >>sys.stderr,"                      a file contains bam filenames and reads number<optional> in each line"
    print >>sys.stderr,"   -n,--normalize     normalize to RPKM"
    print >>sys.stderr,"   -g,--gene          genes annotation table"
    print >>sys.stderr,"   --bed              bed annotation file"
    print >>sys.stderr,"   "

def BamInBed(bamfilename,norm=1000000):
    print >>sys.stderr,"\nProcessing File:",bamfilename
    samfile=pysam.Samfile(bamfilename,"rb")
    reads_number=File_reads_number[bamfilename]
    j=0
    for i in BedList:
        j+=1
        if j%1000==0: print >>sys.stderr,j,"/",len(BedList),"genes","\r",
        try:
            c=samfile.count(i.chr,i.start,i.stop)
        except:
            c=0
        if ifNormalize:
            c=float(c)*norm*1000/reads_number/i.length()  # reads number per 1M 1K bp RPKM?
        if not Table.has_key(i.uid):
            Table[str(i.uid)]=[]
        Table[i.uid].append(c) 
    samfile.close()

def OutputTable():
    
    print "# Bamlist:"
    for i,x in enumerate(bamlist):
        print "#    BAM."+str(i+1)+":",x,"\t",File_reads_number[x],"reads"
    print "# Normalized to RPKM:",ifNormalize
    print "# chr\tstart\tstop\tid\tscore\tstrand",
    for i in range(len(bamlist)):
        print "\tBAM."+str(i+1),
        i
    print 
    for i in BedList:
        print i,
        for j in Table[i.uid]:
            if not ifNormalize:
                print "\t",j,
            else:
                print "\t%.3f"%j,
        print




def Main():
    global Table,Version,ifNormalize,File_reads_number,BedList,bamlist,Norm
    Version="0.1"
    Table={}
    if len(sys.argv)<2:
        show_help()
        exit(0)
    opts,restlist = getopt(sys.argv[1:],"hb:g:n",\
                        ["help","bam=","normalize","gene=","bed="])
    Binsize=2000
    ifNormalize=False
    File_reads_number={}
    ifBamFileList=0
    BedList=[]
    Norm=1000000
    ifUniq=False


    for o, a in opts:                       
        if o in ("-h","--help"): 
            show_help()
            exit
        if o in ("--bam"):
            ifBamFileList=1
            BamFileList=a
        if o in ("--normalize","-n"):
            ifNormalize=True
        if o in ("-g","--gene"):
            genefile=a
        if o in ("--bed"):
            bedfile=a

    if 'genefile' in dir():
        try:
            f=open(genefile)
        except:
            print >>sys.stderr,"can't open genefile",genefile
            exit(0)
        print "# Gene File:",genefile
        for line in f:
            line=line.strip()
            if line[0]=='#': continue
            x=line.split("\t")
            t=GeneBed(x)
            BedList.append(t)
        f.close()
    if 'bedfile' in dir():
        try:
            f=open(bedfile)
            print "# Bed File:",bedfile
            for line in f:
                line=line.strip()
                if line[0]=='#': continue
                x=line.split("\t")
                t=Bed(x)
                BedList.append(t)
            f.close()
        except:
            print >>sys.stderr,"can't open genefile",genefile
            exit(0)
    BedList.sort()
    if len(BedList)==0:
        print >>sys.stderr,"No gene or bed file assigned."
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
        BamInBed(i,Norm)
    OutputTable()

    


    
if __name__=="__main__":
    Main()

