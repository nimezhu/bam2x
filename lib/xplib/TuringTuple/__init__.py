#!/usr/bin/env pythON
# Programmer : zhuxp
# Date: 
# Last-modified: 01-27-2014, 15:11:37 EST
import xplib.Turing.TuringCodeBook as cb
from bitarray import bitarray 
from xplib.Turing.TuringUtils import *
from xplib.Annotation import Bed12


# WARNING: TPDP a lot of things 

# SIMPLE SORTING TUPLES REVISE ABLE

I_POS=0
I_CODE=1
I_CID=2
I_OTHER=3



def turing_tuples_to_graph_string(x,scale=60):
    '''
    translate turing tuples to graph string
    '''
    s=""
    l=x[-1][I_POS]-x[0][I_POS]
    if scale > l:
        scale=l
    for i in range(scale+1):
        s+="-";
    list_s=list(s)
    for i in x:
        position=i[I_POS]-x[0][I_POS]
        index=position*scale/l;
        #print index
        if i[I_CODE]==cb.BLOCKON and list_s[index]=="-": list_s[index]='^'
        if i[I_CODE]==cb.BLOCKOFF and list_s[index]=="-" : list_s[index]='v'
        if i[I_CODE]==cb.BLOCKON and list_s[index]=="v": list_s[index]='x'
        if i[I_CODE]==cb.BLOCKOFF and list_s[index]=="^": list_s[index]='x'
        if i[I_CODE]==cb.ON: list_s[index]="["
        if i[I_CODE]==cb.OFF: list_s[index]=")"
    return "".join(list_s)
 
def codes_length(g):
    h={}
    s=0
    for i in g:
        #print "debug i in g",i
        if i[I_CODE]==cb.BLOCKON or i[I_CODE]==cb.BLOCKOFF:
            if not h.has_key(i[I_POS]):
               h[i[I_POS]]=i[I_CODE]
               s+=1
    return s
def translate_path_into_bits(codes,codes_len,path,bp=5):
    '''
    path and codes should be a sorted list of tuples
    path is the reads
    codes is the transript structure
    small number , so using array instead of TuringTupleSortingArray
    Tuple Sorting List is enough.
    '''
    #TODO: tolerate the reads in +-5bp?
    #TODO  unknown reads if less than given percetage?
    c=list(codes)
    p=list(path)
    c2=list()
    p2=list()
    ccid=0
    pcid=1
    #print "debug",p
    #print "debug",c
    for i in c:
        c2.append(i+(ccid,))
    for i in p:
        #print "debug",i
        # print "debug",pcid
        p2.append(i+(pcid,))
    c=[i for i in heapq.merge(c2,p2)]
    #codes.sort()
    '''
    ccid=c[0].cid
    #print "ccid:",ccid
    pcid=p[0].cid
    #print "pcid:",pcid
    if ccid >= pcid:
        pcid=ccid+1
        for i,x in enumerate(p):
            p[i].cid=pcid
    '''
    cstate=cb.OFF
    cbstate=cb.BLOCKOFF
    pstate=cb.OFF
    pbstate=cb.BLOCKOFF
    #output=bitarray([ True for i in range(2*len(c))])
    output=bitarray(codes_len*2)
    output.setall(True)
    #print "debug",codes
    #print "debug",output


    '''
    revise the mis aligned position for path
    '''
    # UPDATE TUPLE ?
    buff=[]
    last_pos=0
    last_code=cb.ON
    for i,x in enumerate(c):
        if c[i][I_CID]==pcid:
            buff.append(i)
        if c[i][I_CID]==ccid:
            for j in buff:
                if c[j][I_CODE]%2 ==last_code%2 and abs(c[j][I_POS]-last_pos)<=bp:
                    lst=list(c[j])
                    lst[I_POS]=last_pos
                    c[j]=tuple(lst)
                if c[j][I_CODE]%2==c[i][I_CODE]%2 and abs(c[j][I_POS]-c[i][I_POS])<=bp:
                    lst=list(c[j])
                    lst[I_POS]=c[i][I_POS]
                    c[j]=tuple(lst)
            buff=[]
            last_pos=c[i][I_POS]
            last_code=c[i][I_CODE]
    for j in buff:
        if c[j][I_CODE]==last_code and abs(c[j][I_POS]-last_pos)<=bp:
            lst=list(c[j])
            lst[I_POS]=last_pos
            c[j]=tuple(lst)
    '''
    end of revise position
    '''

    j=-1
    last_pos=-1
    # c.sort()
    for i in c:
        #print i[I_CID],i
        if (i[I_POS]!=last_pos):
            if pstate==cb.ON:
                if pbstate==cb.BLOCKON:
                    output[2*j]=True;
                    output[2*j+1]=False;
                if pbstate==cb.BLOCKOFF:
                    if output[2*j]==True and output[2*j+1]==True:
                        output[2*j]=False
            elif pstate==cb.OFF:
                pass
            if(i[I_CID]==ccid and last_pos!=i[I_POS] and (i[I_CODE]==cb.BLOCKON or i[I_CODE]==cb.BLOCKOFF)):
                j+=1
                last_pos=i[I_POS]
        if i[I_CID]==pcid: ## process path
            if i[I_CODE]==cb.ON: 
                pstate=cb.ON
            elif i[I_CODE]==cb.OFF:
                pstate=cb.OFF
                
            elif i[I_CODE]==cb.BLOCKON:
                pbstate=cb.BLOCKON
            elif i[I_CODE]==cb.BLOCKOFF:
                pbstate=cb.BLOCKOFF
        elif i[I_CID]==ccid: ## process self[I_CODE]s
            if i[I_CODE]==cb.ON: 
                cstate=cb.ON
            elif i[I_CODE]==cb.OFF:
                cstate=cb.OFF
            elif i[I_CODE]==cb.BLOCKON:
                cbstate=cb.BLOCKON
            elif i[I_CODE]==cb.BLOCKOFF:
                cbstate=cb.BLOCKOFF
    #print "debug output",output
    #print "return output",output
    return output 
def translate_paths_into_bits(codes,codes_len,paths,bp=5):
    output=bitarray(codes_len*2)
    output.setall(True)
    for i in paths:
        output=bitarray_and(output,translate_path_into_bits(codes,codes_len,i,bp))
    return output

def translate_bits_into_bed(codes,bits,id="noname",chr="chr"):
    '''
    codes is a sorted tuple
    '''
    chr=chr
    id=id
    last_pos=-1
    blockStarts=[]
    blockSizes=[]
    blockCount=0
    itemRgb="0,0,0"
    score=0
    i=0
    exon_state=False
    for x in codes:
        if x[I_POS]!=last_pos:
            last_pos=x[I_POS]
            #print "debug i,",i,bits[2*i],bits[2*i+1]
            #print "debug x",x
            #print "debug last pos:",last_pos
            if bits[2*i]==True and bits[2*i+1]==False:
                if not exon_state and x[I_CODE]==cb.BLOCKON: # TODO : test if and after is right, debut the bug.bed
                    exon_state=True
                    #print "debug adding exon start",last_pos

                    blockStarts.append(last_pos)
                    blockCount+=1
            elif bits[2*i]==False and bits[2*i+1]==True:
                if exon_state:
                    exon_state=False
                    #print "last pos:",last_pos
                    #print "minus:",blockStarts[-1]
                    blockSizes.append(last_pos-blockStarts[-1])
            elif bits[2*i]==True and bits[2*i+1]==True:
                if x[I_CODE]==cb.BLOCKOFF:
                    if exon_state:
                        exon_state=False
                        blockSizes.append(last_pos-blockStarts[-1])
                if x[I_CODE]==cb.BLOCKON:
                    if not exon_state:
                        exon_state=True
                        blockStarts.append(last_pos)
                        blockCount+=1
            i+=1
    if blockCount>0:
        start=blockStarts[0]
        #print "debug",bits
        #print "debug",len(bits)
        #print "debug",codes
        #print "debug",len(codes)
        #print "debug",blockStarts
        #print "debug",blockSizes
        stop=blockStarts[-1]+blockSizes[-1]
        
    else:
        start=0
        stop=0
    cds_start=start
    cds_stop=start
    for i,x in enumerate(blockStarts):
        blockStarts[i]-=start
    strand="+"
    return Bed12([chr,start,stop,id,score,strand,cds_start,cds_stop,itemRgb,blockCount,blockSizes,blockStarts])


##TODO THE TEST


import array,tempfile,heapq
import xplib.Stats.prob as prob
import itertools
from operator import itemgetter
assert array.array('i').itemsize==4

class TuringTupleSortingArray:
    def __init__(self,a=None,MAX_ARRAY_SIZE=500000,**dict):
        self.data=[[]]
        self.files=[]
        self.size=0
        self.index=0
        self.MAX_ARRAY_SIZE=MAX_ARRAY_SIZE
        if a is not None:
            for i in a:
                self.append(i)
        self.sort()
        self.has_sorted=True
        if dict.has_key("tempdir"):
            self.tempdir=dict["tempdir"]
        else:
            self.tempdir=None
    def append(self,x):
        self.has_sorted=False
        if (self.index<self.MAX_ARRAY_SIZE):
            heapq.heappush(self.data[0],x)
            self.index+=1
        else:
            self.data[0].sort(key=itemgetter(0,1))
            f=tempfile.TemporaryFile(dir=self.tempdir)
            self.data.append(TuringTupleSortingArray.tempfile_reader(f))
            TuringTupleSortingArray.write_tempfile(self.data[0],f)
            f.seek(0)
            self.files.append(f)
            self.data[0]=[x]
            self.index=1
    def sort(self):
        self.data[0].sort(key=itemgetter(0,1))
        self.has_sorted=True
    def seek0(self):
        for f in self.files:
            f.seek(0)
    def iter(self):
        #yield "test"
        if not self.has_sorted:
            self.sort()
            self.has_sorted=True
        self.seek0()
        for i in heapq.merge(*self.data):
            #TODO test custom sort?
            yield i
    @staticmethod    
    def tempfile_reader(f):
        while True:
            a = array.array("i")
            a.fromstring(f.read(4000))
            if not a:
                break
            #print str(f),"DEBUG",len(a)
            for i in range(0,len(a),2):
                yield (a[i],a[i+1])
    @staticmethod
    def write_tempfile(a,f):
        b=array.array("i")
        for i in a:
            b.append(i[0])
            b.append(i[1])
        b.tofile(f)
    





def Test():
    a=[(1,cb.ON),(1,cb.BLOCKON),(20,cb.BLOCKOFF),(180,cb.BLOCKON),(200,cb.OFF),(200,cb.BLOCKOFF)]
    b=[(5,cb.ON),(1,cb.BLOCKON),(10,cb.BLOCKOFF),(10,cb.OFF)]
    print turing_tuples_to_graph_string(a)
    print translate_path_into_bits(a,b)
    print translate_path_into_bits(a,a)
        

    
if __name__=="__main__":
    Test()


















