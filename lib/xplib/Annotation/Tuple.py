from collections import namedtuple


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
    
class BED3(namedtuple("BED3",H_BED3),METABED):
    pass


class BED6(namedtuple("BED6",H_BED6),METABED):
    pass


class BED12(namedtuple("BED12",H_BED12),METABED):
    pass
if __name__=="__main__":
    b=("chr1",1,2)
    a=BED3(*b)
    print a
    print a.chr
    print a.end
