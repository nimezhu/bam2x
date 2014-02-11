from collections import namedtuple
import itertools


H_PSL=("matches","misMatches","repMatches","nCount","qNumInsert","qBaseInsert","tNumInsert","tBaseInsert","strand","qName","qSize","qStart","qEnd","tName","tSize","tStart","tEnd","blockCount","blockSizes","qStarts","tStart")
H_BED3=("chr","start","stop")
H_BED6=("chr","start","stop","id","score","strand")
H_BED12=("chr","start","stop","id","score","strand","cds_start","cds_stop","itemRgb","blockCount","blockSizes","blockStarts")
class METABED(object):
    @property
    def end(self):
        return self.stop
    def __cmp__(self,other):
        return cmp(self.chr,other.chr) or cmp(self.start,other.start) or (self.stop,other.stop)
    def __str__(self):
        s=""
        for i in self:
            s+=str(i).strip("(").strip(")")+"\t"
        s.strip("\t")
        s.strip("")
        return s
    def __len__(self):
        return self.stop-self.start
    def cdna_length(self):
        s=0
        for i in self.Exons():
            s=s+len(i)
        return s
    def Exons(self):
        if hasattr(self,"blockCount"):
            a=[]
            if self.strand=="-":
                step=-1
                j=self.blockCount
            else:
                step=1
                j=0
            for i in range(self.blockCount):
                exon_start=self.start+self.blockStarts[i]
                exon_end=self.start+self.blockStarts[i]+self.blockSizes[i]
                j+=step
                exon_id=self.id+"_Exon_"+str(j+1)
                t=(self.chr,exon_start,exon_end,exon_id,self.strand,0.0)
                a.append(BED6(*t))
            if self.strand=="-":
                return a[::-1]
            else:
                return a
        else:
            l=[]
            l.append(self._replace(id=self.id+"_Exon_1"))
            return l
    def Introns(self):
        a=[]
        if hasattr(self,"blockCount"):
            if self.strand=="-":
                step=-1
                j=self.blockCount-1
            else:
                step=1
                j=0
            for i in range(self.blockCount-1):
                intron_start=self.start+self.blockStarts[i]+self.blockSizes[i]
                intron_end=self.start+self.blockStarts[i+1]
                j+=step
                intron_id=self.id+"_Intron_"+str(j+1)
                t=(self.chr,intron_start,intron_end,intron_id,self.strand,0.0)
                a.append(BED6(*t))
            if self.strand=="-":
                return a[::-1]
        return a
    
    
    
class BED3(namedtuple("BED3",H_BED3),METABED):
    pass


class BED6(namedtuple("BED6",H_BED6),METABED):
    pass


class BED12(namedtuple("BED12",H_BED12),METABED):
    def _slice(self,start,end,suffix):
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
        for blockStart,blockSize in itertools.izip(self.blockStarts,self.blockSizes):
            exon_start=blockStart+self.start
            exon_stop=blockStart+blockSize+self.start
            slice_start=max(start,exon_start)
            slice_stop=min(end,exon_stop)
            if slice_start < slice_stop:
                blockCount+=1
                sliceBlockStarts.append(slice_start-start)
                sliceBlockSizes.append(slice_stop-slice_start)
        if blockCount==0:
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




if __name__=="__main__":
    b=("chr1",200,500,"test_gene",0.0,"-",250,450,"0,0,0",2,(100,101),(0,199))
    a=BED12(*b)
    print a
    print a.chr
    print a.end
    for i in a.Exons():
        print i
    for i in a.Introns():
        print i
    print a.utr3()
    print a.utr5()
    print a.cds()
