# Programmer : zhuxp
# Date: 
# Last-modified: 06-30-2014, 15:20:19 EDT

import os,sys
from bam2x.Annotation import *
from bam2x.Annotation import BED6 as Bed
from bam2x.TableIO import hclass
from bam2x import TableIO
from bam2x.Struct import binindex
import pysam
from bam2x.Tools import rc
from bam2x import Tools
from twobitreader import *
import csv
from bam2x import IO
'''
BASIC QUERY FUNCTIONS
'''

def parse_region_str(x):
    try:
        a=x.split(":")
        b=[None,None]
        if len(a)==2:
            b=a.split("-")
            b[0]=int(b[0])-1
            b[1]=int(b[1])
        return BED3(a[0],b[0],b[1])
    except:
        print >>sys.stderr,"unknown formatted region string"
        raise
class MetaDBI:
    '''
    Meta Class for all DBI (Database Interface)
    '''
    def __init__(self,data,**dict):
        self.data=data
    def query(self,x):
        pass
    def __iter__(self):
        for i in self.data:
            yield i
class BinIndexI(MetaDBI):
    '''
    DBI for flat annotation file or list
    cls
    '''
    def __init__(self,file,**dict):
        '''
        Wrapped in bam2x.DBI.init()
        BinIndex(file,cls=inherited_namedtuplecls)
        inherited_namedtuplecls should have _make and _types functino
        or
        assuming the entry in container is already formatted
        BinIndex(container) 
        '''
        if isinstance(file,str):
            file=csv.reader(IO.fopen(file,"r"),delimiter="\t")
        if dict.has_key("cls"):
            cls=dict["cls"]
            if isinstance(cls,str):
                if hclass.has_key(cls):
                    cls=hclass[cls]
                else:
                    print >>sys.stderr,"UNKNOWN FORMAT %s IN BININDEX DATA STRUCT"%cls
    
            self.data=binindex(file,cls=cls)
        else:
            self.data=binindex(file)
    def query(self,x=None,**kwargs):
        '''
        yield the overlap features with x.
        '''
        if x is None:
            try:
                x=BED3(kwargs["chr"],kwargs["start"],kwargs["stop"])
            except:
                print >>sys.stderr,"UNKNOWN QUERY"
                raise 
        if isinstance(x,str):
            x=parse_region_str(x)
        for i in self.data.query(x,**kwargs):
            yield i
    def close(self):
        '''
        release memory
        '''
        self.data=None 
        
class TabixI(MetaDBI):
    '''
    New DBI for tabix file. 
    TabixI(filename,cls=inherited_namedtuplecls)
    '''
    def __init__(self, tabix_file_name,**dict):
        '''
        wrapped in DBI.init(filename,"tabix",cls=cls)
        '''
        self.tabix_file_name=tabix_file_name
        self.dict=dict
        
        if dict.has_key("cls"):
            self.cls=dict["cls"]
            if isinstance(self.cls,str):
                if hclass.has_key(self.cls):
                    self.cls=hclass[self.cls]
                else:
                    print >>sys.stderr,"UNKNOWN FORMAT %s IN BININDEX DATA STRUCT"%cls
        try:
            self.data=pysam.Tabixfile(tabix_file_name)
        except:
            print >>sys.stderr,"WARNING: Can't init the tabix file",tabix_file_name
    def query(self,x=None,**kwargs):
        '''
        yield the overlap feature in tabix index files
        
        dbi=TabixI(filename,cls=cls)
        three options:
        dbi.query(x)  , x is BED class or any class have .chr, .start, .stop
        dbi.query(chr="",start="",stop="")
        dbi.query("chr01:1-10000")
        '''
        try:
            if x is None:
                try:
                    x=BED3(kwargs["chr"],kwargs["start"],kwargs["stop"])
                except:
                    print "UNKNOWN QUERY"
                    raise 
            if isinstance(x,str):
                x=parse_region_str(x)
            if hasattr(self,"cls"):
                for item in self.data.fetch(x.chr,x.start,x.stop):
                    item=self.cls._make(self.cls._types(item.split("\t")))
                    yield item
            else:
                for item in self.data.fetch(x.chr,x.start,x.stop):
                    yield item

        except Exception as e:
           print e
           raise StopIteration
    def close(self):
        self.data.close()
class TwoBitI(MetaDBI): 
    def __init__(self,File,**kwargs):
        '''
        query 2 bit genome
        '''
        self.data=TwoBitFile(File)
    def query(self,x):
        chr=self.data[x.chr]
        return chr[x.start:x.stop]
class GenomeI(TwoBitI):
    '''
    Wrapped for query sequence
    Initialize genome 2bit file only once.
    '''
    def query(self,x=None,**dict):
        if x is None:
            try:
                x=BED3(kwargs["chr"],kwargs["start"],kwargs["stop"])
            except:
                print "UNKNOWN QUERY"
                raise 
        if isinstance(x,str):
            x=parse_region_str(x)
        method="seq"
        if(dict.has_key("method")):
            method=dict["method"]
        if method=="seq":
            return self.get_seq(x)
        elif method=="cdna" or method=="cDNA":
            return self.get_cdna_seq(x)
        elif method=="cds" or method=="CDS":
            return self.get_cds_seq(x)
        elif method=="utr3":
            return self.get_utr3_seq(x)
        elif method=="utr5":
            return self.get_utr5_seq(x)
    def get_seq(self,x):
        chr=self.data[x.chr]
        seq=chr[x.start:x.stop]
        if x.strand=="-":
            seq=rc(seq)
        return seq
    def get_cdna_seq(self,bed12):
        s=""
        for i in bed12.Exons():
            s+=self.get_seq(i)
        return s
    def get_cds_seq(self,bed12):
        cds=bed12.cds()
        if cds is None or len(cds)==0: return ""
        return self.get_cdna_seq(cds)
    def get_utr3_seq(self,bed12):
        utr3=bed12.utr3()
        if utr3 is None or len(utr3)==0: return ""
        return self.get_cdna_seq(utr3)
    def get_utr5_seq(self,bed12):
        utr5=bed12.utr5()
        if utr5 is None or len(utr5)==0: return ""
        return self.get_cdna_seq(utr5)        



class BamlistI(MetaDBI):
    '''
    A DBI for a list of bamfiles ( or a file contain bamfile names) 
    Query a Bed Object (or Annotation Format inherit Bed)
    Yield Nt Distribution from the start position to the end position
    Nt Distribution [A,C,G,T]
    Example:
        [0,12,0,23]  which means  12 C and 23 T 
    Query Example:
        from bam2x import DBI

        bamfiles=["file1.bam","file2.bam","file3.bam"]
        dbi=DBI.init(bamfiles,"bamlist")
        bed=Bed(["chr1",0,5,".",".","."])
        for i,x in enumerate(dbi.query(bed)):
            print "position",i+1,"\t",x

        the output example:
        position 1   [0,1,4,0]
        position 2   [0,0,0,3]
        .
        .
        .
        
        Which indicates that 

        1 C 4 G in position 1
        3 T in position 2 

        ....
        in all file1.bam, file2.bam and file3.bam
        
    '''

    hNtToNum={'a':0,'A':0,
          'c':1,'C':1,
          'g':2,'G':2,
          't':3,'T':3
         }
    Nt=['A','C','G','T']
    
    
    def __init__(self,bamfiles,**dict):
        '''
        '''
        if type(bamfiles)==type("string"):
            filename=bamfiles
            bamfiles=[]
            for i in TableIO.parse(filename,"simple"):
                bamfiles.append(i[0])
        self.bamfiles=[]
        for bamfile in bamfiles:
            if type(bamfile)==type("str"):
                try:
                    bamfile=pysam.Samfile(bamfile,"rb")
                except:
                    print >>sys.stderr,"WARNING: Can't init the bam file",bamfile
            self.bamfiles.append(bamfile)
    def close(self):
        for i in self.bamfiles:
            i.close()
    @property
    def mapped(self):
        s=0
        for i in self.bamfiles:
            s+=i.mapped
        return s
    @property
    def unmapped(self):
        s=0
        for i in self.bamfiles:
            s+=i.unmapped
        return s
    

    def query(self,x=None,method='fetch',**dict):
        try:
            if x is None:
                try:
                    x=BED3(dict["chr"],dict["start"],dict["stop"])
                except:
                    print "UNKNOWN QUERY"
                    raise 
            if type(x)==type("str"):
                x=parse_region_str(x)
            chrom=x.chr
            start=x.start
            end=x.stop
            
            uniq=False
            if dict.has_key("uniq"):
                uniq=dict["uniq"]
            if method=='fetch' or method=="bam2bed12" or method=="fetch12":
                '''
                test version
                still test Tools.cigar_to_coordinates
                '''
                for bamfile in self.bamfiles:
                    for bed in TableIO.parse(bamfile.fetch(chrom,start,end),"bam2bed12",references=chrom,uniq=uniq):
                        yield bed
            elif method=="bam1": 
            # fetch read from paired end with strand information  
                for bamfile in self.bamfiles:
                    strand="read2"
                    if dict.has_key("strand"):   # TODO: if bamfiles have different read1 or read2 ?
                        strand=dict["strand"]
                    for bed in TableIO.parse(bamfile.fetch(chrom,start,end),"bam2bed12",references=chrom,strand=strand,uniq=uniq):
                        yield bed
            elif method=="bam2": # yield bed12
                for bamfile in self.bamfiles:
                    for fragment in TableIO.parse(bamfile.fetch(chrom,start,end),"bam2fragment",bam=bamfile):
                        if dict.has_key("strand"):
                            frag=fragment.toBed12(chr=chrom,strand=dict["strand"],uniq=uniq)
                            if frag: yield frag
                        else:
                            frag=fragment.toBed12(chr=chrom,uniq=uniq)
                            if frag: yield frag
            elif method=="bam2fast": # yield bed12 , don't report mate not in this iterator, 
                for bamfile in self.bamfiles:
                    for fragment in TableIO.parse(bamfile.fetch(chrom,start,end),"bam2fragment"):
                        if dict.has_key("strand"):
                            frag=fragment.toBed12(chr=chrom,strand=dict["strand"],uniq=uniq)
                            if frag: yield frag
                        else:
                            frag=fragment.toBed12(chr=chrom,uniq=uniq)
                            if frag: yield frag
            elif method=='pileup':
                s=[[0,0,0,0] for row in range(end-start)]
                for bamfile in self.bamfiles:
                    try:
                        A=bamfile.pileup(chrom,start,end)
                    except:
                        print >>sys.stderr,"Can't pile up",chrom,start,end
                        raise StopIteration 
                    for pileupcolumn in A:
                        j=pileupcolumn.pos-start
                        if j<0: continue
                        if j>end-start: break
                        for pileupread in pileupcolumn.pileups:
                            try:
                                if pileupread.is_del: continue
                                if pileupread.indel!=0: continue
                                nt=pileupread.alignment.seq[pileupread.qpos]
                                if BamI.hNtToNum.has_key(nt):
                                    k=BamI.hNtToNum[nt]
                                    s[j][k]+=1
                            except:
                                pass
                for i in s:
                    yield i
            elif method=='count':
                s=0
                for bamfile in self.bamfiles:
                    s+=bamfile.count(chrom,start,end)
                yield s
            elif method=='count_fragment':
                s=0
                for bamfile in self.bamfiles:
                   for fragment in TableIO.parse(bamfile.fetch(chrom,start,end),"bam2fragment",bam=bamfile):
                       s+=1
                yield s
            elif method=="seq_and_qual":
                for bamfile in self.bamfiles:
                    for read in bamfile.fetch(chrom,start,end):
                        yield (read.seq,read.qual)
        except:
            raise StopIteration
            

class BamI(BamlistI):
    '''
    A DBI for bamfile
    Query a Bed Object (or Annotation Format inherit Bed)
    Yield Nt Distribution from the start position to the end position

    It is an simple version for BamlistI.
    Use it to  query one BamFile

    Wrapped in DBI.query(bamfilename,"bam")

    Query Example:
        from bam2x import DBI

        bamfiles="file1.bam"
        dbi=DBI.init(bamfiles,"bam")
        bed=Bed(["chr1",0,5,".",".","."])
        for i,x in enumerate(dbi.query(bed)):
            print "position",i+1,"\t",x
    '''
    def __init__(self,bamfile,**dict):
        '''
        init Bam 
        convert filename bamfile into [bamfile].
        '''
        self.bamfiles=[]
        bamfiles=[bamfile]
        for bam in bamfiles:
            if type(bam)==type("str"):
                try:
                    bam=pysam.Samfile(bamfile,"rb")
                except:
                    print >>sys.stderr,"WARNING: Can't init the bam file",bam
            self.bamfiles.append(bam)
 


class BigWigI(MetaDBI):
    '''
    A DBI for bigwig file
    lib dependent: bxpython
    '''
    def __init__(self,bwfile,**dict):
        '''
        init bw file
        '''
        try:
            import bx.bbi.bigwig_file
        except ImportError:
            raise
        if type(bwfile)==type("string"):
            self.data=bx.bbi.bigwig_file.BigWigFile(open(bwfile,"rb"))
        else:
            try:
                self.data=bx.bbi.bigwig_file.BigWigFile(bwfile)
            except:
                print >>sys.stderr,"Error in open bw file"

    def query(self,x=None,**dict):
        '''
        query bw file
        '''
        if x is None:
            try:
                x=BED3(dict["chr"],dict["start"],dict["stop"])
            except:
                print "UNKNOWN QUERY"
                raise 
        if type(x)==type("str"):
            x=parse_region_str(x)
 
        if not dict.has_key("method") or not (dict["method"]=="cDNA" or dict["method"]=="cdna"):
            results=self.data.get_as_array(x.chr,x.start,x.stop)
            if hasattr(x,"strand") and x.strand=="-":
                return results[::-1]
            else:
                return results
        else:
            if dict["method"]=="cDNA" or dict["method"]=="cdna":
                s=[]
                for i in x.Exons():
                    array=self.data.get_as_array(i.chr,i.start,i.stop)
                    size=i.stop-i.start
                    if  x.strand=="+" or x.strand==".":
                        for j in range(size):
                            s.append(array[j])
                    elif x.strand=="-":
                        for j in xrange(size-1,-1,-1):
                            s.append(array[j])
                return s
            else:
                print >>sys.stderr,"query model is not a BED12 class"

        




