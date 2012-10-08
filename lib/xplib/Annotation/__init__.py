# programmer:  zhuxp
# email: nimezhu@163.com
import sys
#Last-modified: 08 Oct 2012 17:53:54
# reader of any column file
__all__=['Bed','GeneBed','TransUnit','Peak','OddsRatioSNP','VCF']        

class Bed(object):
    '''
    Genome Annotation Format.
    using in UCSC genome browser
    Reference : http://genome.ucsc.edu/FAQ/FAQformat.html
    Bed 6
    chrom   start   stop    id  score   strand

    generator example:
    
    bed=Bed(["chr1",1,20,"id",0.0,"+"])
    
    bed=Bed("chr1\t1\t20\tid\t0.0\t+")

    bed=Bed(chr="chr1",start=1,stop=20,id="id",score=0.0,strand="+")
    '''
    def __init__(self,x=None,**kwargs):
        self.id="NONAME"
        self.score=0
        self.strand="."
        if x is not None:
            if type(x)==type("str"):
                x=x.split("\t")
	    self.chr=str(x[0]).strip()
	    self.start=int(x[1])
	    self.stop=int(x[2])
            try:
	        self.id=x[3].strip()
                self.score=float(x[4])
	        self.strand=x[5].strip()
            except:
                pass
        for key in kwargs.keys():
            setattr(self,key,kwargs[key])
    def __str__(self):
        '''
        Return a standard bed 6 format string
        chr start   stop    id  score   strand
        Usage Example:
            a=Bed(chr="chr1",start=1,stop=20)
            print a
        '''
	string=self.chr+"\t"+str(self.start)+"\t"+str(self.stop)+"\t"+str(self.id)+"\t"+str(self.score)
	string+="\t"+str(self.strand)
	return string
    def __len__(self):
        '''
        return the length of the annotation
        '''
        return self.stop-self.start
    def __cmp__(self,other):
	return cmp(self.chr,other.chr) or cmp(self.start,other.start) or cmp(self.stop,other.stop) or cmp(self.strand,other.strand)
    def upstream(self,bp=1000):
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
    def downstream(self,bp=1000):
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
    def tss(self):
        '''return the Bed Object that represent transcription start site, name is geneid_tss'''
        pos=self.stop
        if self.strand=="+":
            pos=self.start
        else:
            pos=self.stop-1  #should minus one or not?
        return Bed([self.chr,pos,pos+1,self.id+"_tss",self.strand,0])
    def tts(self):
        '''return the Bed Object that represent transcription termination site, name is geneid_tss'''
        if self.strand=="+":
            pos= self.stop-1
        else:
            pos=self.start 
        return Bed([self.chr,pos,pos+1,self.id+"_tts",self.strand,0])
    def strand_cmp(self,bed):
        '''
        return if the compare bed is in the same strand or different strand
        "+" means in the same strand
        "-" means in the CR strand
        "." means unknown
        '''
        if bed.strand == "." or self.strand == ".":
            return "."
        if self.strand == bed.strand:
            return "+"
        else:
            return "-"


class GeneBed(Bed):
    '''
    Gene Bed Class
    Download gene annotation table from table browser of UCSC genome Browser
    '''
    def __init__(self,x,**kwargs):
        try:
            self.bin=int(x[0])
            if(self.bin < 10000):
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
        for key in kwargs.keys():
            setattr(self,key,kwargs[key])
        
  #  def __str__(self):
  #      return self.string
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
    def getATG(self):
        if(self.strand == "+"):
            return Bed([self.chr,self.cds_start,self.cds_start+3,self.id+"_ATG",0,self.strand])
        elif (self.strand=="-"):
            return Bed([self.chr,self.cds_stop-3,self.cds_stop,self.id+"_ATG",0,self.strand])
    def getStopCodon(self):
        if(self.strand == "-"):
            return Bed([self.chr,self.cds_start,self.cds_start+3,self.id+"_stop_codon",0,self.strand])
        elif (self.strand=="+"):
            return Bed([self.chr,self.cds_stop-3,self.cds_stop,self.id+"_stop_coden",0,self.strand])
    
    def cds(self):
        '''
        Return the CDS region as a GeneBed Object
        '''
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

    def new_utr5(self):
        '''
        Return 5 UTR as a GeneBed Object
        '''
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
        '''
        Return 3' UTR as a GeneBed Object
        '''
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


                   
class VCF(Bed):
    '''
    VCF OBJECT
    Example:
        for i in TableIO.parse(file,"vcf"):
            print i
    '''
    def __init__(self,x,**kwargs):
        if x is not None:
            sep="\t"
            if kwargs.has_key("sep"): sep=kwargs["sep"]
            if type(x)==type("str"):
                x=x.split(sep)
            self.chr=str(x[0])
            self.chrom=self.chr
            self.pos=int(x[1])
            self.start=self.pos-1
            self.stop=self.pos
            self.id=str(x[2])
            self.ref=str(x[3])
            self.alt=str(x[4])
            try:
                self.qual=float(x[5])
            except:
                self.qual=x[5]
            try:
                self.filter=x[6]
                self.info=x[7]
                self.format=x[8]
                self.others=x[9:]
                self.infos={}
                self._parse_infos()
            except:
                pass
        for key in kwargs.keys():
            setattr(self,key,kwargs[key])

    def _parse_infos(self):
        x=self.info.split(";")
        for i in x:
            a=i.split("=")
            b=a[1].split(",")
            for i,c in enumerate(b):
                try:
                    b[i]=float(c)
                except:
                    b[i]=c
            if len(b)==1:
                self.infos[a[0]]=b[0]
            else:
                self.infos[a[0]]=b
    def getInfo(self,InfoID):
        if self.infos.has_key(InfoID):
            return self.infos[InfoID]
        else:
            return None
    def __str__(self):
        s=""
        s+=str(self.chr)+"\t"+str(self.pos)+"\t"+self.id+"\t"+self.ref+"\t"+self.alt
        s+="\t"+str(self.qual)
        if self.filter is not None:
            s+="\t"+str(self.filter)
        if self.info is not None:
            s+="\t"+str(self.info)
        if self.format is not None:
            s+="\t"+str(self.format)
        try:
            for i in self.others:
                s+="\t"+str(i)
        except:
            pass
        return s
        
    def __cmp__(self,other):
	return cmp(self.chr,other.chr) or cmp(self.start,other.start)


class Repeat(Bed):
    def __init__(self,x=None,**kwargs):
        sep="\t"
        if kwargs.has_key("sep"):
            sep=kwargs["sep"]
        if type(x)==type(str):
            x=x.split(sep)
        if x is not None:
            try:
                self.bin=int(x[0])
                self.swScore=int(x[1])
                self.milliDiv=int(x[2])
                self.milliDel=int(x[3])
                self.millIns=int(x[4])
                self.chr=str(x[5])
                self.genoName=self.chr
                self.start=int(x[6])
                self.genoStart=self.start
                self.stop=int(x[7])
                self.genoEnd=self.stop
                self.genoLeft=int(x[8])
                self.strand=str(x[9])
                self.repName=str(x[10])
                self.repClass=str(x[11])
                self.repFamily=str(x[12])
                self.repStart=int(x[13])
                self.repEnd=int(x[14])
                self.repLeft=int(x[15])
                self.id=x[16]
            except:
                print >>sys.stderr,"wrong format file"
        for key in kwargs.keys():
            setattr(self,key,kwargs[key])
    def __str__(self):
        s=""
        s+=str(self.bin)
        s+="\t"+str(self.swScore)
        s+="\t"+str(self.milliDiv)
        s+="\t"+str(self.milliDel)
        s+="\t"+str(self.millIns)
        s+="\t"+str(self.genoName)
        s+="\t"+str(self.genoStart)
        s+="\t"+str(self.genoEnd)
        s+="\t"+str(self.genoLeft)
        s+="\t"+str(self.strand)
        s+="\t"+str(self.repName)
        s+="\t"+str(self.repClass)
        s+="\t"+str(self.repFamily)
        s+="\t"+str(self.repStart)
        s+="\t"+str(self.repEnd)
        s+="\t"+str(self.repLeft)
        s+="\t"+str(self.id)
        return s






###################### Below is Private Format
class Peak(Bed):
    '''
    ChIP Seq Peak Class
    extension Bed Class
    Private Class
    read the output of xbams2peak.py
    '''
    def __init__(self,x,**kwargs):
        self.chr=x[0].strip()
        self.chr=self.chr.replace("'","")
        self.start=int(x[1])
        self.end=int(x[2])
        self.stop=int(x[2])   
        self.reads_num=int(x[3])
        self.pvalue=float(x[4])
        self.coverage=float(x[5])
        self.peak_pos=int(x[6])
        self.peak_coverage=int(x[7])
    def __str__(self):
        s="("
        s+=self.chr+","
        s+=str(self.start)+","
        s+=str(self.end)+","
        s+=str(self.reads_num)+","
        s+=str(self.pvalue)+","
        s+=str(self.coverage)+","
        s+=str(self.peak_pos)+","
        s+=str(self.peak_coverage)+")"
        return s
    def __cmp__(self,other):
	return cmp(self.chr,other.chr) or cmp(self.start,other.start) or cmp(self.stop,other.stop) 
    def tab(self):
        s=self.chr+"\t"
        s+=str(self.start)+"\t"
        s+=str(self.end)+"\t"
        s+=str(self.reads_num)+"\t"
        s+=str(self.pvalue)+"\t"
        s+=str(self.coverage)+"\t"
        s+=str(self.peak_pos)+"\t"
        s+=str(self.peak_coverage)
        return s
class TransUnit(Bed):
    '''
    Private File Format
    For parse the output of xbams2trans.py

    Usage:
        for i in TableIO.parse("File.Trans","trans"):
            print i
    Or  
        for i in TableIO.parse("File.Trans","trans"):
            a.append(i)
        a.sort()
        for i in a:
            print i
    '''
    def __init__(self,x=None):
        self.genes=[]
        self.feats=[]
        self.promoters=[]
        self.promoterInfo=""
        if x:
            if type(x)==type("s"):
                x=x.split("\t")
            self.processHeaer(x)
    def processHeader(self,x):
        if x[0]=="NV " or x[0]=="NP " or x[0]=="KN " or x[0]=="NV" or x[0]=="NP" or x[0]=="KN":
            self.group=x[0]
            x=x[1:]
        self.name=x[0]
        self.chr=x[1]
        self.start=int(x[2])
        self.end=int(x[3])
        self.stop=int(x[3])   
        self.reads_num=int(x[4])
        self.pvalue=float(x[5])
        self.coverage=float(x[6])
        self.peak_pos=int(x[7])
        self.peak_coverage=int(x[8])
        self.guess_strand=x[9]
    def __str__(self):
        s=""
        s+=self.header()+"\n"
        s+=self.group+"\t"+"Promoter Info:\t"+self.promoterInfo+"\n"
        for g in self.genes:
            s+=self._str_gene(g)+"\n"
        for f in self.feats:
            s+=self._str_feat(f)+"\n"
        for p in self.promoters:
            s+=self._str_promoter(p)+"\n"
        s+="//\n"
        return s
    def _str_gene(self,g):
        return self.group+"\t"+"OverlapGene:\t"+str(g)
    def _str_feat(self,f):
        return self.group+"\tOverlapFeat:\t"+str(f)
    def _str_promoter(self,p):
        return self.group+"\tNearbyPromoter:\t"+str(p)

    def __cmp__(self,other):
	return cmp(self.chr,other.chr) or cmp(self.start,other.start) or cmp(self.stop,other.stop) 

    def header(self):
        s=""
        s+=self.group+"\t"
        s+=self.name+"\t"
        s+=self.chr+"\t"
        s+=str(self.start)+"\t"
        s+=str(self.end)+"\t"
        s+=str(self.reads_num)+"\t"
        s+=str(self.pvalue)+"\t"
        s+=str(self.coverage)+"\t"
        s+=str(self.peak_pos)+"\t"
        s+=str(self.peak_coverage)+"\t"
        s+=self.guess_strand
        return s
    def append_overlap_gene(self,gene):
        self.genes.append(gene)
    def append_overlap_feat(self,feat):
        self.feats.append(feat)
    def append_promoter(self,promoter):
        if type(promoter)==type("s"):
            s,strand=promoter.split(")")
            s=s.replace("(","")
            s=s.replace(",","\t")
            x=s.split("\t")
            self.promoters.append(Peak(x))
        elif type(promoter)==type([1,2,3]):
            self.promoters.append(Peak(promoter))
        
class OddsRatioSNP(Bed):
    '''
    OddsRatio Table are generated by program xbams2APS.py
    position are 0-index 
    Example Format:
    chr1    10086   A/G     3.75495232253   ( 21 1 14 6 )   [20, 1, 0, 2] vs [13, 0, 5, 1]
    
    chromosome  
    position    
    MajorAllele/MinorAllele 
    APS Score       
    (Major in Case,Minor in Case,Major in Control, Minor in Control)
    [A,C,G,T] in Case vs [A,C,G,T] in Control
    '''
    def __init__(self,x=None,**kwargs):
        '''
        initialize object 
        
        Examples:
            x=OddsRatioSNP([chr,start,stop,major_allele,.....])

            x=OddsRatioSNP(chr="chr1",start=pos-1,A=[2,0,4,0],B=[4,0,1,0])
        APS Score :
            Allele Preference Score
                Chi2 Value For Log Odds Ratio
            Reference:
                http://en.wikipedia.org/wiki/Odds_ratio
        '''
        if x is None:
            self.chr=None
            self.start=None
            self.stop=None
            self.major_allele=None
            self.minor_allele=None
            self.odds_ratio=None
            self.APS=None
            self.odds_ratio_matrix=None
            self.A_nt_dis=None
            self.B_nt_dis=None
        else:
            if type(x)==type("str"):
                x=x.split("\t")
            self.chr=x[0].strip()
            self.start=int(x[1]) # 0-index for OddsRatio Table
            self.stop=int(x[1])+1
            a=x[2].split("/")
            self.major_allele=a[0].strip()
            self.minor_allele=a[1].strip()
            self.odds_ratio=float(x[3]) #old format compatible
            self.APS=float(x[3])
            x[4]=x[4].replace("( ","")
            x[4]=x[4].replace(" )","")
            self.odds_ratio_matrix=x[4].split(" ")
            for i in range(4):
                self.odds_ratio_matrix[i]=int(self.odds_ratio_matrix[i])
            a=x[5].split("vs")
            for i in range(2):
                a[i]=eval(a[i])
            self.A_nt_dis=a[0]
            self.B_nt_dis=a[1]
        binA=False
        binB=False
        if kwargs.has_key("A"):
            self.A_nt_dis=kwargs['A']
            binA=True
        if kwargs.has_key('B'):
            self.B_nt_dis=kwargs['B']
            binB=True
        if binA and binB:
            self.calculate_from_nt_dis()
        if kwargs.has_key("chr"):
            self.chr=kwargs["chr"]
        if kwargs.has_key("start"):
            self.start=kwargs['start']
            self.stop=self.start+1
        

    def __str__(self):
        '''
        return the standard format 
        Usage:
            print aps
        '''
        s=""
        s+=self.chr+"\t"+str(self.start)+"\t"+self.major_allele+"/"+self.minor_allele+"\t"+str(self.APS)
        s+="\t"+"( "
        for i in self.odds_ratio_matrix:
            s+=str(i)+" "
        s+=")\t"+str(self.A_nt_dis)+" vs "+str(self.B_nt_dis)
        return s
        
    def __cmp__(self,other):
	return cmp(self.chr,other.chr) or cmp(self.start,other.start) or cmp(self.stop,other.stop) 
    def calculate_from_nt_dis(self):
        '''
        calculate the aps score from nt distribution
        '''
        from math import log
        s=[0,0,0,0]
        idx=[0,0,0,0]
        Nt=['A','C','G','T']
        a0=self.A_nt_dis
        a1=self.B_nt_dis
        for i in range(4):
            s[i]=a0[i]+a1[i]
        for i in range(4):
            for j in range(i+1,4):
                if s[i] < s[j]: idx[i]+=1
                if s[i] >= s[j]: idx[j]+=1
        for i in range(4):
            if idx[i]==0: idx1=i
            if idx[i]==1: idx2=i
        (a11,a12,a21,a22)=(a0[idx1],a0[idx2],a1[idx1],a1[idx2])  # zero to one
        if a11==0 or a12==0 or a21==0 or a22==0:
            ratio=(float(a11+0.5)/float(a12+0.5))/(float(a21+0.5)/float(a22+0.5))
            logratio=log(ratio)
            sigma2=1.0/(a11+1)+1.0/(a12+1)+1.0/(a21+1)+1.0/(a22+1)
        else:
            ratio=(float(a11)/float(a12))/(float(a21)/float(a22))
            logratio=log(ratio)
            sigma2=1.0/a11+1.0/a12+1.0/a21+1.0/a22
        x=logratio*logratio/sigma2
        self.APS=x
        self.major_allele=Nt[idx1]
        self.minor_allele=Nt[idx2]
        self.odds_ratio_matrix=(a11,a12,a21,a22)
        if logratio > 0: self.A_enrich=self.major_allele
        elif logratio < 0: self.A_enrich=self.minor_allele
        else: self.A_enrich=None
