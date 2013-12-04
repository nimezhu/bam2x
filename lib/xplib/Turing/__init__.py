#!/usr/bin/env pythON
# Programmer : zhuxp
# Date: 
# Last-modified: 12-04-2013, 17:21:30 EST
import TuringCodeBook as cb
from bitarray import bitarray 
from TuringUtils import *
from xplib.Annotation import Bed12
class TuringCode:
    def __init__(self,pos,code,**kwargs):
        self.pos=pos
        self.code=code
        self.cid=cb.DEFAULT_CID
        for key in kwargs.keys():
            setattr(self,key,kwargs[key])
    def __str__(self):
        s=str(self.pos)+":"+str(self.code)
        return s
    def __cmp__(self,other):
        return cmp(self.pos,other.pos) or cmp(self.code,other.code)




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

    def iter_all_path(self,i):
        #print self.states_register 
        if self.codes[i].code==cb.ON:
            self.state=cb.BLOCKOFF
            self.states_register.append(i)
            for x in self.iter_all_path(i+1): yield x
        elif self.state==cb.BLOCKOFF and self.codes[i].code==cb.OFF:  ## OFF < BLOCKOFF
   #     elif self.codes[i].code==cb.OFF:  ## OFF < BLOCKOFF
            x = list(self.states_register)
            x.append(i)
            yield x
        elif self.codes[i].code==cb.BLOCKON and self.state==cb.BLOCKOFF:
            self.states_register.append(i)
            self.state=cb.BLOCKON
            for x in self.iter_all_path(i+1): yield x
            self.states_register.pop()
            self.state=cb.BLOCKOFF
            for x in self.iter_all_path(i+1): yield x
        elif self.codes[i].code==cb.BLOCKOFF and self.state==cb.BLOCKON:
            self.states_register.append(i)
            self.state=cb.BLOCKOFF
            for x in self.iter_all_path(i+1): yield x
            self.states_register.pop()
            self.state=cb.BLOCKON
            for x in self.iter_all_path(i+1): yield x
        elif self.codes[i].code==cb.BLOCKON and self.state==cb.BLOCKON:
            for x in self.iter_all_path(i+1): yield x
        elif self.codes[i].code==cb.BLOCKOFF and self.state==cb.BLOCKOFF:
            for x in self.iter_all_path(i+1): yield x
    def iter_turing_paths(self):
        for i in self.iter_all_path(0):
            l=[]
            for j in i:
                l.append(TuringCode(self.codes[j].pos,self.codes[j].code))
            g2=TuringGraph(l)
            yield g2
    
    def count_paths_number(self):
        '''
        test version
        count the number of paths
        '''
        count_on=0
        count_off=1
        for code in self.codes:
            if code.code==cb.BLOCKON:
               count_on=count_on+count_off 
            elif code.code==cb.BLOCKOFF:
               count_off=count_off+count_on
        return count_off 


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
        d={}
        for i in self.codes:
            if d.has_key(i.pos):
                d[i.pos]+=1
            else:
                d[i.pos]=1
        return len(d.values())
    def score_path(path):
        #TO DO
        pass
    def translate_path_into_bits(self,path):
        #TO DO
        c=list(self.codes)
        p=list(path.codes)
        ccid=c[0].cid
        #print "ccid:",ccid
        pcid=p[0].cid
        #print "pcid:",pcid
        if ccid==pcid:
            pcid+=1
            for i,x in enumerate(p):
                p[i].cid=pcid
        cstate=cb.OFF
        cbstate=cb.BLOCKOFF
        pstate=cb.OFF
        pbstate=cb.BLOCKOFF
        #output=bitarray([ True for i in range(2*len(c))])
        output=bitarray(len(self)*2)
        output.setall(True)
        
        for i in p:
            c.append(i)
        c.sort()
        j=-1
        last_pos=-1

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
    def translate_paths_into_bits(self,paths):
        a=[]
        output=bitarray(len(self)*2)
        output.setall(True)
        for i in paths:
            output=bitarray_and(output,self.translate_path_into_bits(i))
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
        return Bed12([chr,start,stop,id,score,strand,cds_start,cds_stop,itemRgb,blockCount,blockSizes,blockStarts])







def test():
    a=list()
    codes=list()
    h={}
    codes.append( TuringCode(150,cb.BLOCKON))
    codes.append( TuringCode(1,cb.ON))
    codes.append( TuringCode(1,cb.BLOCKON))
    codes.append( TuringCode(100,cb.BLOCKOFF))
    codes.append( TuringCode(200,cb.BLOCKON))
    codes.append( TuringCode(300,cb.OFF))
    codes.append( TuringCode(300,cb.BLOCKOFF))
    codes.sort()

    paths=list()
    paths.append( TuringCode(1,cb.ON))
    paths.append( TuringCode(1,cb.BLOCKON))
    paths.append( TuringCode(100,cb.BLOCKOFF))
    paths.append( TuringCode(200,cb.BLOCKON))
    paths.append( TuringCode(300,cb.BLOCKOFF))
    paths.append( TuringCode(300,cb.OFF))
    g= TuringGraph(codes)
    p= TuringGraph(paths)

    print g;
    print g.graph_str();
    for i in g.iter_all_path(0):
        print "yield:",i
        l=[]
        for j in i:
            print j
            l.append(TuringCode(g.codes[j].pos,g.codes[j].code))
        g2=TuringGraph(l)
        print g2
        print g2.graph_str()
        print g.translate_path_into_bits(g2)
        print g.translate_bits_into_bed(g.translate_path_into_bits(g2))
    #    a.append(g.translate_path_into_bits(TuringGraph(l)))
    #    h[g.translate_path_into_bits(TuringGraph(l))]=1

    for x in g.codes: print x
    print g.translate_path_into_bits(p)

    print g.translate_bits_into_bed(g.translate_path_into_bits(p))

    '''
    for i,x in enumerate(a):
        for j,y in enumerate(a):
            print i,j
            print x,y
            print twobitarray_and(x,y)
    print h
    '''




if __name__=="__main__":
    test()



