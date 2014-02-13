# Programmer : zhuxp
# Date: 
# Last-modified: 02-12-2014, 22:02:57 EST
import os,sys,argparse,types
from bam2x.Annotation import BED6 as Bed 
from bam2x import TableIO
'''
binindex Data Structure

self.data is hash table. 
Key is Chromosome
self.data["chr1"] is a list of feature list
self.data["chr"][bin_number] is a list of features that have bin index.

The Most Useful Function is:

Read features into a data structure:
    Example:
        f=TableIO.parse("filename","bed12")
        data=binindex(f)
    Example:
        f=TableIO.parse("filename","vcf")
        data=binindex(f)
    Example:
        data=binindex(bedlist)
    

Query the overlap features in the data structure:
    Usage:
        for i in data.query(bed):
            print i

This BinIndex Data Structure was wrapped into DBI.init() and DBI.query():
    Usage Example:
        dbi=DBI.init("file.bed","bed")
        x=Bed(chr="chr1",start=1,stop=200)
        for i in dbi.query(x):
            print i
TO DO:
    Now BinIndex only have item inherent from Bed Class, which means it need to have at least three attributes, chr, start and stop. and start is 0-index.  
    it use the fourth attribute, id.
    Need to improve to read any table that contain this information but not not in a pre defined format. 
'''
class binindex(object):
	'''
        binindex data structure
        this data structure hold genome annotations. 
        it can be queried overlap features with given chromosome regions.
        
        Examples:
        
        Query Overlap Features:
        
        data=binindex(TableIO.parse("file.bed","bed"))
        bed=Bed(chr="chr1",start=1,stop=10000)
        for i in data.query(bed):
            print i
        Add Two Data Struture:
            dataC = dataA + dataB
        
        '''
	binOffsets=(512+64+8+1,64+8+1,8+1,1,0)
	binFirstShift=17
	binNextShift=3
        binLength=(4096+512+64+8+1)
        @staticmethod
	def range2bin(start,end):
	    '''
	        Convert (start,end) to bin index
	        
	        Usage:
	           data=[[] for row in range(4096+512+64+8+1)]
	           for i in list:
	               bin=range2bin(i.start,i.stop)
	               data[bin].append(i)
	    '''
	    startBin=start
	    endBin=end-1 
	    startBin >>= binindex.binFirstShift
	    endBin >>= binindex.binFirstShift
	    for i in range(len(binindex.binOffsets)):
	        if startBin==endBin:
	            return binindex.binOffsets[i]+startBin
	        startBin >>= binindex.binNextShift
	        endBin >>= binindex.binNextShift
	    print >>sys.stderr,"start %d , end %d out of range"%(start,end) 
	    return None
            
        @staticmethod
	def iter_range_overlap_bins(start,end):
	    '''
	        Iterate the possible overlap binindex
	        the core function of query overlap features.
	        Usage:
	            for i in binindex.iter_range_overlap_bins(start,end):
	                for j in data[i]:
	                    if start > j.stop and stop < j.start:
	                        print j
	    '''
	    startBin=start
	    endBin=end-1
	    startBin >>= binindex.binFirstShift
	    endBin >>= binindex.binFirstShift
	    for i in range(len(binindex.binOffsets)):
	        for j in range(startBin,endBin+1):
	            yield j+binindex.binOffsets[i]
	        startBin >>= binindex.binNextShift
	        endBin >>= binindex.binNextShift
        @staticmethod
        def bin2cov(data):
            '''
            input a binindex data
            output total length of coverage for each bin 
            '''
            coverage={}
            for index,i in enumerate(data):
                if index%100==0:
                    print >>sys.stderr,"processing %d entries\r"%index,
                if not coverage.has_key(i.chr):
                    coverage[i.chr]=[0 for k in range(binindex.binLength)]
                i_bin=binindex.range2bin(i.start,i.stop)
                for j in binindex.iter_range_overlap_bins(i.start,i.stop):
                    (start,stop)=binindex.bin2range(j)
                    max_start=max(start,i.start)
                    min_stop=min(stop,i.stop)
                    length=min_stop-max_start
                    if float(length)/(stop-start)<0:
                        print >>sys.stderr,length,stop,start
                    coverage[i.chr][j]+=length
            return coverage
        
        
        @staticmethod
        def bin2range(bin):
            binShift=binindex.binFirstShift
            for i in binindex.binOffsets:
                if bin-i>=0:
                    bin=bin-i
                    break
                binShift+=binindex.binNextShift
            return (bin<<binShift,(bin+1)<<binShift)
        
        @staticmethod
        def bin2length(bin):
            (start,end)=binindex.bin2range(bin)
            return end-start
        @staticmethod
        def bin2level(bin):
            for (i,x) in enumerate(binindex.binOffsets):
                if bin-x>=0:
                    return 4-i 

      
        
        def __init__(self,handle=None,**kwargs):
            self.data={}
            if handle is not None:
                self.read(handle,**kwargs)
            
	def append(self,bed):
	    '''
	    data is a binindex Data Structure
	    bed is an anntation feature. ( bed vcf or genebed etc.)
	    append an annotation into data
	    '''
	    a=bed
	    bin=binindex.range2bin(a.start,a.stop)
	    if not self.data.has_key(a.chr):
	        self.data[a.chr]=[[] for row in range(binindex.binLength)]
	    self.data[a.chr][bin].append(a)
	
	def remove(self,bed):
	    '''
	    delete bed entry from a BinIndex Data Structure data.
	    '''
	    a=bed
	    bin=binindex.range2bin(a.start,a.stop)
	    for i,x in enumerate(self.data[a.chr][bin]):
	        if a==x:                          # use define equal in Bed.__cmp__
	            del self.data[a.chr][bin][i]
	
	def merge(self,other):
	    '''
	    merge two binindex data structure
	    and only ratain one if duplicates.
	    '''
	    for i in other:
	        flag=1
	        for j in self.query(i):
	            if j.start==i.start and j.stop==i.stop and j.id==i.id:
	                flag=0
	        if flag:
	            self.append(i)
        def __add__(self,b):
            '''
            overload + 
            merge two data structure
            and retain all information including duplicate entries.
            '''
            c=binindex()
            for i in self:
                c.append(i)
            for i in b:
                c.append(i)
            return c

	
	    
	def __iter__(self):
	    '''
	    iter all the features in binindex structure
	    Usase:
	        data is a binindex data structure
	
	        for i in iterBinIndex(data):
	            print i
            should revise to iter in sorted order
	    '''
	    chrs=self.data.keys()
	    chrs.sort()
	    for chr in chrs:
	        for binindex in range(4096+512+64+8+1):
	            for i in self.data[chr][binindex]:
	                yield i
	
	def read(self,handle,**kwargs):
	    '''
	    Read iterator into BinIndex Data Structure
	    if have cls, turn tuple into namedtuple
            '''
            cls=None
            if kwargs.has_key("cls"):
                cls=kwargs["cls"]
	    if cls is not None:
                for i in handle:      
                    a=cls._make(cls._types(i))
	            self.append(a)
            else:
                for i in handle:
                    self.append(i)
	
	def query(self,bed,**kwargs):
	    '''
	    iterator the bed overlap features in data
	    Usage:
                data=bedindex(file.bed)
	        for i in data.query(bed):
	            print i
	
	    **kwargs for further extension for annotation with no pre define format
	    '''
	
	    if type(bed)==type((1,2,3)) or type(bed)==([1,2,3]):
	        bed=Bed(*bed[0:3])  # guess (chrome,start,stop)
	    if not self.data.has_key(bed.chr):
	         raise StopIteration 
	    D=self.data[bed.chr]
	    for bin in binindex.iter_range_overlap_bins(bed.start,bed.stop):
	        for f in D[bin]:
	            if f.start < bed.stop and f.stop > bed.start:
	                yield f
	    raise StopIteration
        def __len__(self):
            s=0
            for chr in self.data:
                for j in self.data[chr]:
                    s+=len(j)
            return s
        def __str__(self):
            s="binindex data structure object with "+str(len(self))+" entries"
            return s
        def uniq(self):
            '''
            return a binindex data structure that have no duplicate entries.

            Usage:
            dataUniq=data.uniq()
            '''
            c=binindex()
            for i in self:
                if c.contains(i): continue
                else:
                    c.append(i)
            return c
        def contains(self,a):
            '''
            binary decision that if this data struct has feature a.
            '''
            if not self.data.has_key(a.chr):
                return False
	    bin=binindex.range2bin(a.start,a.stop)
	    for i,x in enumerate(self.data[a.chr][bin]):
	        if a==x:
                    return True
            return False


	
	
	
def test():
    if len(sys.argv)==1:
        print >>sys.stderr,"Usage: binindex.py file.bed"
        exit()
    a=TableIO.parse(sys.argv[1],'bed12')
    data=binindex(a)
    data2=binindex()
    bed=Bed("chr1",100000,2000000,".",".",".")
    for i in data.query(bed):
        print "before remove:",len(data)
        data.remove(i)
        print "after remove:",len(data)
        data2.append(i)
        print data2
    for i in data2:
        print i
    print "data finalize:"
    data.merge(data2)
    print data
    print data+data2
    print data
    print data.uniq()
    print data

    
    
if __name__=="__main__":
    test()


