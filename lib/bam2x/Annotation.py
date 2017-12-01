from collections import namedtuple
import itertools
import logging

H_PSL=("matches","misMatches","repMatches","nCount","qNumInsert","qBaseInsert","tNumInsert","tBaseInsert","strand","qName","qSize","qStart","qEnd","tName","tSize","tStart","tEnd","blockCount","blockSizes","qStarts","tStart")
H_VCF=("chr","pos","id","ref","alt","qual","filter","info")
F_VCF=(str,int,str,str,str,int,str,str)
H_BED3=("chr","start","stop")
F_BED3=(str,int,int)
H_BED4=("chr","start","stop","id")
F_BED4=(str,int,int,str)
H_BED6=("chr","start","stop","id","score","strand")
F_BED6=(str,int,int,str,float,str)
H_BED12=("chr","start","stop","id","score","strand","cds_start","cds_stop","itemRgb","blockCount","blockSizes","blockStarts")
def int_tuple(x):
    a=[]
    for i in x.strip().strip(",").split(","):
        a.append(int(i))
    return tuple(a)
F_BED12=(str,int,int,str,float,str,int,int,str,int,int_tuple,int_tuple)
class METABED(object):
    @property
    def end(self):
        return self.stop
    def __cmp__(self,other):
        return cmp(self.chr,other.chr) or cmp(self.start,other.start) or cmp(self.stop,other.stop)
    def __str__(self):
        return "\t".join([str(i) for i in self])
    def __len__(self):
        return self.stop-self.start
    def Exons(self):
        if hasattr(self,"blockCount"):
            a=[]
            if self.strand=="-":
                step=-1
                j=self.blockCount-1
            else:
                step=1
                j=0
            for i in range(self.blockCount):
                exon_start=self.start+self.blockStarts[i]
                exon_end=self.start+self.blockStarts[i]+self.blockSizes[i]
                exon_id=self.id+"_Exon_"+str(j+1)
                j+=step
                t=(self.chr,exon_start,exon_end,exon_id,0.0,self.strand)
                a.append(BED6(*t))
            if self.strand=="-":
                return a[::-1]
            else:
                return a
        else:
            l=[]
            if hasattr(self,"id"):
                l.append(self._replace(id=self.id+"_Exon_1"))
            else:
                l.append(self)
            return l
    def Introns(self):
        a=[]
        if hasattr(self,"blockCount"):
            if self.strand=="-":
                step=-1
                j=self.blockCount-2
            else:
                step=1
                j=0
            for i in range(self.blockCount-1):
                intron_start=self.start+self.blockStarts[i]+self.blockSizes[i]
                intron_end=self.start+self.blockStarts[i+1]
                intron_id=self.id+"_Intron_"+str(j+1)
                j+=step
                t=(self.chr,intron_start,intron_end,intron_id,0.0,self.strand)
                a.append(BED6(*t))
            if self.strand=="-":
                return a[::-1]
        return a
    def cdna_length(self):
        '''sum of the length of exons'''
        s=0
        for i in self.Exons():
            s=s+i.stop-i.start
        return s

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
        if (start<0):
            start=0
            id=self.id+"up"+str(self.start)
        x=[chr,start,stop,id,0,strand]
        return BED6(*x)
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
        return BED6(*x)
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
        if (start<0):
            start=0
            id=self.id+"down"+str(self.start)
        x=[chr,start,stop,id,0,strand]
        return BED6(*x)
    def tss(self):
        '''return the Bed Object that represent transcription start site, name is geneid_tss'''
        pos=self.stop
        if self.strand=="+":
            pos=self.start
        else:
            pos=self.stop-1  #should minus one or not?
        return BED6(self.chr,pos,pos+1,self.id+"_tss",0,self.strand)
    def tts(self):
        '''return the Bed Object that represent transcription termination site, name is geneid_tss'''
        if self.strand=="+":
            pos= self.stop-1
        else:
            pos=self.start
        return BED6(self.chr,pos,pos+1,self.id+"_tts",0,self.strand)
    def head(self,bp=2):
        '''return the head 2(default) bp of simple bed'''
        if self.strand=="+":
            return BED6(self.chr,self.start,self.start+bp,self.id+"_head"+str(bp)+"bp",0,self.strand)
        if self.strand=="-":
            return BED6(self.chr,self.stop-bp,self.stop,self.id+"_head"+str(bp)+"bp",0,self.strand)
    def tail(self,bp=2):
        '''return the tail 2(default) bp of simple bed'''
        if self.strand=="+":
            return BED6(self.chr,self.stop-bp,self.stop,self.id+"_tail"+str(bp)+"bp",0,self.strand)
        if self.strand=="-":
            return BED6(self.chr,self.start,self.start+bp,self.id+"_tail"+str(bp)+"bp",0,self.strand)
def types(x,TYPES):
    b=[]
    for v,t in itertools.izip(x,TYPES):
        b.append(t(v))
    return tuple(b)
class BED3(namedtuple("BED3",H_BED3),METABED):
    @classmethod
    def _types(cls,x):
        return types(x,F_BED3)
    @property
    def strand(self):
        return "."
    @property
    def id(self):
        return "noname"
    @property
    def score(self):
        return 0.0
class BED4(namedtuple("BED4",H_BED4),METABED):
    @classmethod
    def _types(cls,x):
        return types(x,F_BED4)
    @property
    def strand(self):
        return "."
    @property
    def score(self):
        return 0.0
class BED6(namedtuple("BED6",H_BED6),METABED):
    @classmethod
    def _types(cls,x):
        return types(x,F_BED6)
    @property
    def cds_start(self):
        return self.start
    @property
    def cds_stop(self):
        return self.stop
    @property
    def iterRgb(self):
        return "0,0,0"
    @property
    def blockCount(self):
        return 1
    @property
    def blockStarts(self):
        return [0]
    @property
    def blockSizes(self):
        return [self.stop-self.start]
class VCF(namedtuple("VCF",H_VCF),METABED):
    @classmethod
    def _types(cls,x):
        return types(x,F_VCF)
    @property
    def start(self):
        return self.pos-1
    @property
    def stop(self):
        return self.pos


class BED12(namedtuple("BED12",H_BED12),METABED):
    def __cmp__(self,other):
        return cmp(self.chr,other.chr) or cmp(self.start,other.start) or cmp(self.stop,other.stop) or cmp(self.strand,other.strand) or cmp(self.blockCount,other.blockCount) or cmp(",".join([str(i) for i in self.blockSizes]),",".join([str(i) for i in other.blockSizes])) or cmp(",".join([str(i) for i in self.blockStarts]),",".join([str(i) for i in other.blockStarts]))
    def _slice(self,start,end,suffix="sliced"):
        chr=self.chr

        if start < self.start: start=self.start
        if end > self.stop: end=self.stop
        strand=self.strand
        id=self.id+"_"+suffix
        score=self.score
        itemRgb=self.itemRgb
        cds_start=max(start,self.cds_start)
        cds_stop=min(end,self.cds_stop)
        blockCount=0
        sliceBlockStarts=[]
        sliceBlockSizes=[]

        sliceStart = end
        sliceEnd = start
        for blockStart,blockSize in itertools.izip(self.blockStarts,self.blockSizes):
            exon_start=blockStart+self.start
            exon_stop=blockStart+blockSize+self.start
            slice_start=max(start,exon_start)
            slice_stop=min(end,exon_stop)
            if slice_start < slice_stop:
                if sliceStart > slice_start:
                    sliceStart = slice_start
                if sliceEnd < slice_stop:
                    sliceEnd = slice_stop
                blockCount+=1
                sliceBlockStarts.append(slice_start-start)
                sliceBlockSizes.append(slice_stop-slice_start)
        if start < sliceStart:
            offset = sliceStart - start
            start = sliceStart
            for i,v in enumerate(sliceBlockStarts):
                sliceBlockStarts[i] -= offset
        if end > sliceEnd:
            offset = end - sliceEnd
            end = sliceEnd
        if blockCount==0:
            #logging.warn("wrong slice {start} to {end} for gene {g}".format(start=start,end=end,g=self))
            return None
        else:
            return BED12(chr,start,end,id,score,strand,cds_start,cds_stop,itemRgb,blockCount,tuple(sliceBlockSizes),tuple(sliceBlockStarts))
    def utr5(self):
        if self.strand=="+":
            return self._slice(self.start,self.cds_start,"utr5")
        elif self.strand=="-":
            return self._slice(self.cds_stop,self.stop,"utr5")
        else:
            return None
    def utr3(self):
        if self.strand=="-":
            return self._slice(self.start,self.cds_start,"utr3")
        elif self.strand=="+":
            return self._slice(self.cds_stop,self.stop,"utr3")
        else:
            return None
    def cds(self):
        return self._slice(self.cds_start,self.cds_stop,"cds")
    def __str__(self):
        s=""
        for i,x in enumerate(self):
            if i==10 or i==11:
                for y in x:
                    s+=str(y)+","
                s+="\t"
            else:
                s+=str(x)+"\t"
        return s.strip("\t")
    @classmethod
    def _types(cls,x):
        return types(x,F_BED12)
    @property
    def exon_starts(self):
        return [int(i+self.start) for i in self.blockStarts]

    @property
    def exon_stops(self):
        return [int(i+j+self.start) for i,j in itertools.izip(self.blockStarts,self.blockSizes)]


class Fragment:
    '''
    Paired End Raw Data
    '''
    def __init__(self,read,mate=None,**kwargs):
        self.reads=[]
        if mate is None:
            self.reads.append(read)
        else:
            if read.is_read1:
                self.reads.append(read)
                self.reads.append(mate)
            else:
                self.reads.append(mate)
                self.reads.append(read)
        if kwargs.has_key("chr"):
            self.chr=kwargs["chr"]
    def __str__(self):
        s=""
        for i,x in enumerate(self.reads):
            s+="READ"+str(i)+"\t"+str(x)+"\n"
        return s
    def toBed12(self,chr="unknown_chr",strand="read2",**dict):
        from bam2x import TableIO
        x=list()
        for i in TableIO.parse(self.reads,"bam2bed12",references=chr,strand=strand,**dict):
            x.append(i)
        return x


if __name__=="__main__":
    from bam2x.Tools import translate
    b=("chr1",200,500,"test_gene",0.0,"-",250,450,"0,0,0",2,(100,101),(0,199))
    a=BED12(*b)
    print a
    print a.chr
    print a.end
    for i in a.Exons():
        print i
        for j in translate(a,i):
            print j
    for i in a.Introns():
        print i
        for j in translate(a,i):
            print j
    print a.utr3()
    print a.utr5()
    print a.cds()
    print a.exon_starts[1]
    print a.exon_stops[1]
