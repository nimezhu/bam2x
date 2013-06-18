# Programmer : zhuxp
# Date: 
# Last-modified: 06-18-2013, 16:43:33 EDT

import os,sys
from xplib.Annotation import *
from xplib import TableIO
from xplib.Struct import binindex
import pysam
from xplib.Tools import rc
from twobitreader import *
'''
BASIC QUERY FUNCTIONS
'''
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
    Read annotations (bed,genebed or vcf etc.) into binindex structure.
    Please read documents in xplib.Annotation.Utils for BinIndex Structure detail
    Query the overlap features using xplib.Annotation.Utils.iterOverlapFeature.
    Usage:
            from xplib import DBI
            dbi=DBI.init(file handle or filename or list,"bed")
            for i in dbi.query(bed):
                print i
    '''
    def __init__(self,file,**dict):
        '''
        Wrapped in xplib.DBI.init()
        '''
        if type(file)==type([1,2,3]):
            f=file
        else:
            format=dict['format']
            f=TableIO.parse(file,format)
        self.data=binindex(f)
    def query(self,x):
        '''
        yield the overlap features with x.
        '''
        for i in self.data.query(x):
            yield i
        
class TabixI(MetaDBI):
    '''
    DBI for tabix file.
    it is a simple interface for pysam.Tabix.
    Usage:
        from xplib import DBI
        dbi=DBI.init(filename,"tabix")
        for i in dbi.query(bed):
            print i

    '''
    def __init__(self, tabix_file_name,**dict):
        '''
        wrapped in DBI.init(filename,"tabix")
        '''
        self.tabix_file_name=tabix_file_name
        try:
            self.data=pysam.Tabixfile(tabix_file_name)
        except:
            print >>sys.stderr,"WARNING: Can't init the tabix file",tabix_file_name
        if dict.has_key("header") and dict["header"]==True:
            f=TableIO.parse(tabix_file_name)
            h=f.next()
            l=len(h)
            for i in range(l):
                h[i]=h[i].strip()
            self.header=h
            f.close()


    def query(self,x):
        '''
        yield the overlap feature in tabix index files
        '''
        try:
            for item in self.data.fetch(x.chr,x.start,x.stop):
                yield item
        except:
           raise StopIteration
class FormatTabixI():
    def __init__(self,format,tabixfile):
        pass
    def query(self,x):
        pass
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
    def get_seq(self,bed):
        seq=self.query(bed)
        if bed.strand=="-":
            seq=rc(seq)
        return seq
    def get_cdna_seq(self,genebed):
        s=""
        for i in genebed.Exons():
            s+=self.get_seq(i)
        return s
    def get_cds_seq(self,genebed):
        cds=genebed.cds()
        if cds is None or len(cds)==0: return ""
        return self.get_cdna_seq(cds)
    def get_utr3_seq(self,genebed):
        utr3=genebed.utr3()
        if utr3 is None or len(utr3)==0: return ""
        return self.get_cdna_seq(utr3)
    def get_utr5_seq(self,genebed):
        utr5=genebed.utr5()
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
        from xplib import DBI

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

    def query(self,x,method='pileup'):
        if method=='fetch':
            for bamfile in self.bamfiles:
                for read in bamfile.fetch(x.chr,x.start,x.stop):
                    if read.tid<0:continue
                    strand='+'
                    if read.is_reverse:
                        strand='-'
                    score=read.mapq
                    bed=Bed([bamfile.references[read.tid],read.pos,read.aend,read.qname,score,strand])
                    yield bed
        elif method=='pileup':
            s=[[0,0,0,0] for row in range(x.stop-x.start)]
            for bamfile in self.bamfiles:
                try:
                    A=bamfile.pileup(x.chr,x.start,x.stop)

                except:
                    print >>sys.stderr,"Can't pile up",x.chr,x.start,x.stop
                    raise StopIteration 
                for pileupcolumn in A:
                    j=pileupcolumn.pos-x.start
                    if j<0: continue
                    if j>x.stop-x.start: break
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



class BamI(BamlistI):
    '''
    A DBI for bamfile
    Query a Bed Object (or Annotation Format inherit Bed)
    Yield Nt Distribution from the start position to the end position

    It is an simple version for BamlistI.
    Use it to  query one BamFile

    Wrapped in DBI.query(bamfilename,"bam")

    Query Example:
        from xplib import DBI

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
 


import bx.bbi.bigwig_file
class BigWigI(MetaDBI):
    '''
    A DBI for bigwig file
    lib dependent: bxpython
    '''
    def __init__(self,bwfile,**dict):
        '''
        init bw file
        '''
        if type(bwfile)==type("string"):
            self.data=bx.bbi.bigwig_file.BigWigFile(open(bwfile,"rb"))
        else:
            try:
                self.data=bx.bbi.bigwig_file.BigWigFile(bwfile)
            except:
                print >>sys.stderr,"Error in open bw file"

    def query(self,x,**dict):
        '''
        query bw file
        '''
        if not dict.has_key("model"):
            return self.data.get_as_array(x.chr,x.start,x.stop)
        else:
            if dict["model"]=="cDNA" or dict["model"]=="cdna":
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
                print >>sys.stderr,"query model is wrong for bw file"

        




