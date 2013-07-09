#!/usr/bin/python
# Zhu Xiaopeng, <nimezhu@gmail.com>; Wei Guifeng, <guifengwei@gmail.com> some modification
import sys,re,string
from math import log,sqrt
from zSeqIO import *
# Last-modified: 13 Jun 2013 12:48:34 PM

class GeneBed(Bed):
    ''' GeneClass: GenePred-format
        add the .promoter(bp=1000)
        the exon and intron will be adjusted if you want to exlude the region overlaped with promoter
            Bedlist    self.adjust_exons
            Bedlist    self.adjust_introns
            Bed        self.adjust_utr3
        You can also access the exons,introns and utr through .exons(), .introns(), .utr5(), .utr3()
    '''
    def __init__(self,x):
        try:
            self.bin=int(x[0])
            if(self.bin < 10000):
                x=x[1:]
        except:
            pass
        self.id=x[0]
        self.chr=x[1]
        self.strand=x[2]
        if(self.strand == 1):
            self.strand="+"
        elif(self.strand == -1):
            self.strand="-"
        elif(self.strand == 0):
            self.strand="+"
        self.start=int(x[3])
        self.stop=int(x[4])
        self.cds_start=int(x[5])
        self.cds_stop=int(x[6])
        self.exon_count=int(x[7])
        self.exon_starts=x[8].split(",")
        self.exon_stops=x[9].split(",")
        for i in range(self.exon_count):
            try:
                self.exon_starts[i]=int(self.exon_starts[i])
                self.exon_stops[i]=int(self.exon_stops[i])
            except:
                pass
        self.score=0
        try:
            self.protein_id=x[11]
            self.name2=x[11]
        except:
            self.name2="None"
        try:
            self.align_id=x[12]
        except:
            pass
        self.exons=self.Exons()
        self.introns=self.Introns()
    def __str__(self):
        return "%s\t%s\t%c\t%d\t%d\t%d\t%d\t%d\t%s,\t%s,\t%s\t%s" % (self.id,self.chr,self.strand,self.start,self.stop,self.cds_start,self.cds_stop,self.exon_count,",".join([str(self.exon_starts[j]) for j in range(self.exon_count)]),",".join([str(self.exon_stops[j]) for j in range(self.exon_count)]),str(self.score),self.name2)
    def _exon(self,i):
        '''internal fucntion to call the exon position'''
        if i > self.exon_count:
            return None
        if self.strand=="+" :
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
    def get_cdna_seq(self,fn="/home/wyf/ws190/ce6.2bit"):
        s=""
        for i in self.Exons():
            s+=i.getSeq(fn)
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
    def utr5(self):   #need to change to UTR bed class
        if(self.strand == "+"):
            if(self.cds_start == self.start):
                return None
            return Bed([self.chr,self.start,self.cds_start,self.id+"_"+"utr5",0,self.strand])
        if(self.strand == "-"):
            if(self.cds_stop == self.stop):
                return None
            return Bed([self.chr,self.cds_stop,self.stop,self.id+"_"+"utr5",0,self.strand]) #need to change to UTRBed instead of Bed
    def utr3(self):  #need to change to UTR bed class
        if(self.strand == "-"):
            if(self.cds_start == self.start):
                return None
            return Bed([self.chr,self.start,self.cds_start,self.id+"_"+"utr3",0,self.strand])
        if(self.strand == "+"):
            if(self.cds_stop == self.stop):
                return None
            return Bed([self.chr,self.cds_stop,self.stop,self.id+"_"+"utr3",0,self.strand])

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
    def promoter(self,bp=1000):
        ''' Return the promoter region of the gene(Bed format) '''
        if not bp:
            bp = self.bp
        self.adjust_structure()
        if self.strand == '+':
            return Bed([self.chr,self.start-bp, self.start+bp, 'promoter_'+str(bp)+'_of_'+str(self.id),0,'+'])
        else:
            return Bed([self.chr,self.stop-bp, self.stop+bp, 'promoter_'+str(bp)+'_of_'+str(self.id),0,'-'])
    def adjust_structure(self):
        ''' Depend on the promoter region, the exons,introns,utr3 are adjusted '''
        exons=self.Exons()
        introns=self.Introns()
        utr3=self.utr3()
        a,b,c = [],[],[]
        if self.strand == '+':
            adjust_tss = self.start + self.bp
            for i in exons:
                if i.start < adjust_tss < i.stop or i.start >= adjust_tss:
                    start = max(i.start, adjust_tss)
                    a.append( Bed([i.chr, start, i.stop, 'adj_'+i.id,0,'+']) )
            for i in introns:
                if i.start < adjust_tss < i.stop:
                    b.append( Bed( [i.chr, adjust_tss, i.stop, 'adj_'+i.id,0,'+']) )
                if i.start >= adjust_tss:
                    b.append( Bed( [i.chr, i.start, i.stop, 'adj_'+i.id,0,'+']) )
            if utr3:
                if adjust_tss < utr3.stop:
                    c=Bed( [utr3.chr, max(utr3.start, adjust_tss), utr3.stop, 'adj_'+utr3.id, 0, '+'] )
        else:
            adjust_tss = self.stop - self.bp
            for i in exons:
                if i.stop <= adjust_tss or i.stop > adjust_tss >= i.start:
                    stop = min(i.stop, adjust_tss)
                    a.append( Bed([i.chr, i.start, stop, 'adj_'+i.id,0,'-']) )
            for i in introns:
                if i.stop <= adjust_tss or i.stop > adjust_tss >= i.start:
                    stop = min(i.stop, adjust_tss)
                    b.append( Bed([i.chr, i.start, stop, 'adj_'+i.id,0,'-']) )
            if utr3:
                if adjust_tss > utr3.start:
                    c=Bed([utr3.chr, utr3.start, min(adjust_tss,utr3.stop),'adj_'+utr3.id, 0,'-'])
        self.adjust_exons = a
        self.adjust_introns = b
        self.adjust_utr3 = c
