#!/usr/bin/python
# nimezhu@163.com
import sys
#Last-modified: 29 Jun 2012 23:46:56

# reader of any column file
        

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
        return Bed([self.chr,pos,pos+1,self.id+"_tss",self.strand,0])
    def tts(self):
        if self.strand=="+":
            pos= self.stop-1
        else:
            pos=self.start 
        return Bed([self.chr,pos,pos+1,self.id+"_tts",self.strand,0])
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
        # self.string="\t".join(x)
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
    def plant_promoter(self):
        if(self.strand == "+"):
            return Bed([self.chr,self.cds_start,self.cds_start-2000,self.id+"_ATG",0,self.strand])
        elif (self.strand=="-"):
            return Bed([self.chr,self.cds_stop,self.cds_stop+2000,self.id+"_ATG",0,self.strand])
    def getATGStart(self):
        if(self.strand == "+"):
            return self.cds_start
        elif(self.strand == "-"):
            return self.cds_stop
    def getStopCodonPos(self):
        if(self.strand == "-"):
            return self.cds_start
        elif(self.strand == "+"):
            return self.cds_stop
    def distance_to_atg(self,bed):
        if bed.chr != self.chr : 
            return None
        atg=self.getATGStart()
        pos=(bed.start+bed.stop)/2
       # strand=self.strand_cmp(bed)
        if(self.strand=="+"):
            rpos=pos-atg
        else:
            rpos=atg-pos
        return rpos
    def distance_to_stop_codon(self,bed):
        if bed.chr != self.chr : 
            return None
        scp=self.getStopCodonPos()
        pos=(bed.start+bed.stop)/2
       # strand=self.strand_cmp(bed)
        if(self.strand=="+"):
            rpos=pos-scp
        else:
            rpos=scp-pos
        return rpos
    def bed_in_gene_body_percent(self,bed):
        if bed.chr != self.chr :
            return None
        pos=(bed.start+bed.stop)/2
        length=self.length()
        if(self.strand=="+"):
            return (pos-self.start)*100/length
        else:
            return (self.stop-pos)*100/length
    def bed_in_CDS_percent(self,bed):
        if bed.chr != self.chr :
            return None
        pos=(bed.start+bed.stop)/2
        length=self.cds_stop-self.cds_start
        if(self.strand=="+"):
            return (pos-self.cds_start)*100/length
        else:
            return (self.cds_stop-pos)*100/length
    def bed_gene_report_V1(self,bed):
        if(self.distance_to_atg(bed)<0):
            return "UP",self.distance_to_atg(bed)
        elif(self.distance_to_stop_codon(bed)>0):
            return "DOWN",self.distance_to_stop_codon(bed)
        else:
            return "BODY",self.bed_in_CDS_percent(bed)
        
    def PrintString(self):
        print self.chr,
        print self.id,self.exon_starts,self.exon_stops
    def cds(self):
        return self.new_cds()
    def utr5(self):
        return self.new_utr5()
    def utr3(self):
        return self.new_utr3()
    
    def new_cds(self): 
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
    def new_utr3(self): 
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
    def bed_mid_in_cDNA_percent(self,bed):
        if(self.chr!=bed.chr): return -1
        pos=(bed.start+bed.stop)/2
        cdna_len=self.cdna_length()
        s=0
        for e in self.Exons():
            if self.strand=="+":
                if e.stop > pos:
                    s=s+pos-e.start
                    break
                else:
                    s=s+e.stop-e.start
            if self.strand=="-":
                if e.start < pos:
                    s=s+e.stop-pos
                    break
                else:
                    s=s+e.stop-e.start
        if s<0: return None
        if s>cdna_len: return None
        return float(s)/cdna_len

                    
        
