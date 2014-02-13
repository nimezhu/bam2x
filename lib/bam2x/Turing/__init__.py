#!/usr/bin/env pythON
# Programmer : zhuxp
# Date: 
# Last-modified: 02-13-2014, 10:03:34 EST
import TuringCodeBook as cb
from bitarray import bitarray 
from TuringUtils import *
from bam2x.Annotation import BED12 as Bed12
__all__=["TuringCode","TuringGraph","TuringSortingArray"]

H_TURING3=("pos","code","cid")
H_TURING2=("pos","code")

from collections import namedtuple


class TuringCode3(namedtuple("TuringCode3",H_TURING3)):
    def __str__(self):
        s=str(self.pos)+":"+str(self.code)+ " cid "+str(self.cid)
        return s
    def __cmp__(self,other):
        return cmp(self.pos,other.pos) or cmp(self.cid,other.cid) or cmp(self.code,other.code)


class TuringCode2(namedtuple("TuringCode2",H_TURING2)):
    def __str__(self):
        s=str(self.pos)+":"+str(self.code)
        return s
    def __cmp__(self,other):
        return cmp(self.pos,other.pos) or cmp(self.code,other.code)

TuringCode=TuringCode2

class TuringGraph:
    def __init__(self,codes,paths=[],**kwargs):
        self.codes=codes
        self.paths=paths
        self.state=cb.BLOCKON
        self.states_register=[]
        self.cid=cb.DEFAULT_CID
        self.codes.sort()
        for key in kwargs.keys():
            setattr(self,key,kwargs[key])
    def set_id(self,id):
        self.id=id
    def get_id(self):
        if self.__dict__.has_key("id"):
            return self.id
        else:
            return "noname"

    def set_cid(self,cid):
        self.cid=cid
        for i in range(len(self.codes)):
            self.codes[i].cid=cid



    def add(path):
        #TO DO
        self.paths.append(path)
        pass
    def __str__(self):
        s=""
        for i in self.codes:
            s+=str(i.pos)+":"+str(i.code)+"\n"
        return s
    def graph_str(self,scale=60):
        s=""
        l=self.codes[-1].pos-self.codes[0].pos
        if scale > l:
            scale=l
        for i in range(scale+1):
            s+="-";
        list_s=list(s)
        
        for i in self.codes:
            position=i.pos-self.codes[0].pos
            index=position*scale/l;
            #print index
            if i.code==cb.BLOCKON and list_s[index]=="-": list_s[index]='^'
            if i.code==cb.BLOCKOFF and list_s[index]=="-" : list_s[index]='v'
            if i.code==cb.BLOCKON and list_s[index]=="v": list_s[index]='x'
            if i.code==cb.BLOCKOFF and list_s[index]=="^": list_s[index]='x'
        return "".join(list_s)
    def path_str(self,path,scale=60):
        s=""
        l=self.codes[-1].pos-self.codes[0].pos
        if scale > l:
            scale=l
        for i in range(scale+1):
            s+=" ";
        list_s=list(s)
        
        for i in path.codes:
            position=i.pos-self.codes[0].pos
            index=position*scale/l;
            print index
            if i.code==cb.BLOCKON and list_s[index]==" ": list_s[index]='^'
            if i.code==cb.BLOCKOFF and list_s[index]==" " : list_s[index]='v'
            if i.code==cb.BLOCKON and list_s[index]==" ": list_s[index]='x'
            if i.code==cb.BLOCKOFF and list_s[index]==" ": list_s[index]='x'
        return "".join(list_s)
 

    def __len__(self):
        if hasattr(self,"_len"):
            return self._len
        d={}
        for i in self.codes:
            if d.has_key(i.pos):
                d[i.pos]+=1
            else:
                d[i.pos]=1
        self._len=len(d.values())
        return self._len
    def score_path(path):
        #TO DO
        pass
    def translate_path_into_bits(self,path,bp=5):
        #TODO: tolerate the reads in +-5bp?
        #TODO  unknown reads if less than given percetage?
        c1=list(self.codes)
        p1=list(path.codes)
        ccid=1
        pcid=2
        c=list()
        for i in c1:
            c.append(TuringCode3(i.pos,i.code,ccid))
        for i in p1:
            c.append(TuringCode3(i.pos,i.code,pcid))
        c.sort()
        cstate=cb.OFF
        cbstate=cb.BLOCKOFF
        pstate=cb.OFF
        pbstate=cb.BLOCKOFF
        #output=bitarray([ True for i in range(2*len(c))])
        output=bitarray(len(self)*2)
        output.setall(True)

        '''
        revise the mis aligned position for path
        '''
        buff=[]
        last_pos=0
        last_code=cb.ON
        for i,x in enumerate(c):
            if c[i].cid==pcid:
                buff.append(i)
            if c[i].cid==ccid:
                for j in buff:
                    if c[j].code%2 ==last_code%2 and abs(c[j].pos-last_pos)<=bp:
                        c[j]=c[j]._replace(pos=last_pos)
                    if c[j].code%2==c[i].code%2 and abs(c[j].pos-c[i].pos)<=bp:
                        c[j]=c[j]._replace(pos=c[i].pos)
                buff=[]
                last_pos=c[i].pos
                last_code=c[i].code
        for j in buff:
            if c[j].code==last_code and abs(c[j].pos-last_pos)<=bp:
                c[j]=c[j]._replace(pos=last_pos)
        '''
        end of revise position
        '''
        j=-1
        last_pos=-1
        #c.sort()
        for i in c:
            #print i.cid,i
            if (i.pos!=last_pos):
                if pstate==cb.ON:
                    if pbstate==cb.BLOCKON:
                        output[2*j]=True;
                        output[2*j+1]=False;
                    if pbstate==cb.BLOCKOFF:
         #               print pbstate,"BLOCKOFF"
                        if output[2*j]==True and output[2*j+1]==True:
                            output[2*j]=False
                elif pstate==cb.OFF:
                    pass
                if(i.cid==ccid):
                    j+=1
                    last_pos=i.pos
            if i.cid==pcid: ## process path
                if i.code==cb.ON: 
                    pstate=cb.ON
                elif i.code==cb.OFF:
                    pstate=cb.OFF
                    
                elif i.code==cb.BLOCKON:
                    pbstate=cb.BLOCKON
                elif i.code==cb.BLOCKOFF:
                    pbstate=cb.BLOCKOFF
            elif i.cid==ccid: ## process self.codes
                if i.code==cb.ON: 
                    cstate=cb.ON
                elif i.code==cb.OFF:
                    cstate=cb.OFF
                elif i.code==cb.BLOCKON:
                    cbstate=cb.BLOCKON
                elif i.code==cb.BLOCKOFF:
                    cbstate=cb.BLOCKOFF
        return output 
    def translate_paths_into_bits(self,paths,bp=5):
        a=[]
        output=bitarray(len(self)*2)
        output.setall(True)
        for i in paths:
            output=bitarray_and(output,self.translate_path_into_bits(i,bp))
        return output
    def translate_bits_into_bed(self,bits,id="noname"):
        chr=self.get_id()
        last_pos=-1
        blockStarts=[]
        blockSizes=[]
        blockCount=0
        itemRgb="0,0,0"
        score=0
        id=id
        i=0
        exon_state=False
        for x in self.codes:
            if x.pos!=last_pos:
                last_pos=x.pos
                #print "debug last pos:",last_pos
                if bits[2*i]==True and bits[2*i+1]==False:
                    if not exon_state:
                        exon_state=True
                        blockStarts.append(last_pos)
                        blockCount+=1
                elif bits[2*i]==False and bits[2*i+1]==True:
                    if exon_state:
                        exon_state=False
                        #print "last pos:",last_pos
                        #print "minus:",blockStarts[-1]
                        blockSizes.append(last_pos-blockStarts[-1])
                elif bits[2*i]==True and bits[2*i+1]==True:
                    if x.code==cb.BLOCKOFF:
                        if exon_state:
                            exon_state=False
                            blockSizes.append(last_pos-blockStarts[-1])
                    if x.code==cb.BLOCKON:
                        if not exon_state:
                            exon_state=True
                            blockStarts.append(last_pos)
                            blockCount+=1
                i+=1
        if blockCount>0:
            start=blockStarts[0]
            stop=blockStarts[-1]+blockSizes[-1]
        else:
            start=0
            stop=0
        cds_start=start
        cds_stop=start
        for i,x in enumerate(blockStarts):
            blockStarts[i]-=start
        strand="+"
        return Bed12(chr,start,stop,id,score,strand,cds_start,cds_stop,itemRgb,blockCount,blockSizes,blockStarts)


    def translate_bits_into_beds(self,bits,id="noname",chr="chr"):
        '''
        beds means to fragmentation the bits.
        not try to link them with missing data
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
        i0=0
        strand="+"
        for x in self.codes:
            if x.code!=cb.BLOCKON and x.code!=cb.BLOCKOFF: continue
            if x.pos!=last_pos:
                last_pos=x.pos
                if bits[2*i]==True and bits[2*i+1]==False:
                    if not exon_state and x.code==cb.BLOCKON: # TODO : test if and after is right, debut the bug.bed
                        exon_state=True
                        blockStarts.append(last_pos)
                        blockCount+=1
                elif bits[2*i]==False and bits[2*i+1]==True:
                    if exon_state:
                        exon_state=False
                        blockSizes.append(last_pos-blockStarts[-1])
                elif bits[2*i]==True and bits[2*i+1]==True:
                    if blockCount>0:
                        start=blockStarts[0]
                        if exon_state:
                            exon_state=False
                            blockSizes.append(last_pos-blockStarts[-1])
                        stop=blockStarts[-1]+blockSizes[-1]
                        for j,y in enumerate(blockStarts):
                            blockStarts[j]-=start
                        i0+=1
                        yield Bed12(chr,start,stop,id+"."+str(i0),score,strand,start,start,itemRgb,blockCount,blockSizes,blockStarts)
                        blockCount=0
                        blockStarts=[]
                        blockSizes=[]
                i+=1
    
        if blockCount>0:
            if exon_state:
                blockSizes.append(last_pos-blockStarts[-1])
                exon_state=False
            start=blockStarts[0]
            stop=blockStarts[-1]+blockSizes[-1]
            for j,y in enumerate(blockStarts):
                blockStarts[j]-=start
            i0+=1
            yield Bed12(chr,start,stop,id+"."+str(i0),score,strand,start,start,itemRgb,blockCount,blockSizes,blockStarts)
                    

import array,tempfile,heapq
import itertools
from operator import itemgetter
assert array.array('i').itemsize==4

class TuringSortingArray:
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
            self.data[0].sort()
            f=tempfile.TemporaryFile(dir=self.tempdir)
            self.data.append(TuringSortingArray.tempfile_reader(f))
            TuringSortingArray.write_tempfile(self.data[0],f)
            f.seek(0)
            self.files.append(f)
            self.data[0]=[x]
            self.index=1
    def sort(self):
        self.data[0].sort()
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
                yield TuringCode2(a[i],a[i+1])
    @staticmethod
    def write_tempfile(a,f):
        b=array.array("i")
        for i in a:
            b.append(i.pos)
            b.append(i.code)
        b.tofile(f)

    




def test():
    a=list()
    codes=list()
    h={}
    codes.append( TuringCode2(150,cb.BLOCKON))
    codes.append( TuringCode2(1,cb.ON))
    codes.append( TuringCode2(1,cb.BLOCKON))
    codes.append( TuringCode2(100,cb.BLOCKOFF))
    codes.append( TuringCode2(200,cb.BLOCKON))
    codes.append( TuringCode2(300,cb.OFF))
    codes.append( TuringCode2(300,cb.BLOCKOFF))
    codes.sort()

    paths=list()
    paths.append( TuringCode2(1,cb.ON))
    paths.append( TuringCode2(1,cb.BLOCKON))
    paths.append( TuringCode2(100,cb.BLOCKOFF))
    paths.append( TuringCode2(200,cb.BLOCKON))
    paths.append( TuringCode2(300,cb.BLOCKOFF))
    paths.append( TuringCode2(300,cb.OFF))
    g= TuringGraph(codes)
    p= TuringGraph(paths)

    print g;
    print g.graph_str();

    for x in g.codes: print x
    print g.translate_path_into_bits(p)
    print bitarray_to_rep(g.translate_path_into_bits(p))

    print g.translate_bits_into_bed(g.translate_path_into_bits(p))





if __name__=="__main__":
    test()



