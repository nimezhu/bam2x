from xplib.Annotation import Bed12
import itertools
tuple_header={"chr":0,"start":1,"stop":2,"id":3,"score":4,"strand":5}
CHROM=0
CHROMSTART=1
CHROMEND=2
NAME=3
SCORE=4
STRAND=5
THICKSTART=6
THICKEND=7
ITEMRGB=8
BLOCKCOUNT=9
BLOCKSIZES=10
BLOCKSTARTS=11
def tuple_cmp(a,b):
    '''
    cmp two bed tuples
    '''
    return cmp(a[CHROM],b[CHROME]) or cmp(a[CHROMSTART],b[CHROMSTART]) or cmp(a[CHROMEND],b[CHROMEND])
def tuple_factory(a):
    return Bed12(a)
def tuple_sort(l):
    l.sort(cmp=tuple_cmp)
def tuple_len(a):
    return a[CHROMEND]-a[CHROMSTART]
'''
below is function for bed
'''
def get_exons(a):
    '''
    return bed tuple or a list of bed tuple?
    or an iterator?
    '''
    #TODO
    l=[]
    if a[STRAND]=="-":
        for i,x in enumerate(zip(reversed(a[BLOCKSTARTS]),reversed(a[BLOCKSIZES]))):
            l.append((a[CHROM],a[CHROMSTART]+x[0],a[CHROMSTART]+x[0]+x[1],a[NAME]+"_Exon_"+str(i),a[STRAND],0.0))
    else:
        for i,x in enumerate(zip(a[BLOCKSTARTS],a[BLOCKSIZES])):
            l.append((a[CHROM],a[CHROMSTART]+x[0],a[CHROMSTART]+x[0]+x[1],a[NAME]+"_Exon_"+str(i),a[STRAND],0.0))
    
    return l
def get_introns(a):
    '''
    same as get_exons
    to test!
    '''
    #TODO
    l=[]
    if a[STRAND]=="+":
        for i in range(a[BLOCKCOUNT]-1):
            l.append((a[CHROM],a[CHROMSTART]+a[BLOCKSTARTS][i]+a[BLOCKSIZES][i],a[CHROMSTART]+a[BLOCKSTARTS][i+1],a[NAME]+"_Intron_"+str(i),a[STRAND],0.0))
    elif a[STRAND]=="-":
        for i in range(a[BLOCKCOUNT]-1):
            l.append((a[CHROM],a[CHROMSTART]+a[BLOCKSTARTS][i]+a[BLOCKSIZES][i],a[CHROMSTART]+a[BLOCKSTARTS][i+1],a[NAME]+"_Intron_"+str(a[BLOCKCOUNT]-2-i),a[STRAND],0.0))
    return l

def is_overlap(a,b):
    '''
    '''
    if a[CHROM]==b[CHROM] and a[START] < b[STOP] and b[START] < a[STOP] :
        return True
    else:
        return False
def get_upstream(a):
    '''
    '''
    #TODO 
    pass
def get_downstream(a):
    '''
    '''
    #TODO
    pass

def test():
    print "test"
    bed=("chr01",1,323,"test1",0.0,"-",3,320,"0,0,0",2,(30,272),(0,50))
    print bed
    for i in get_exons(bed):
        print i
    for i in get_introns(bed):
        print i


if __name__=="__main__":
    test()
