#!/usr/bin/env pythON
# Programmer : zhuxp
# Date: 
# Last-modified: 11-12-2013, 16:04:00 EST
import TuringCodeBook as cb
from bitarray import bitarray 
from TuringUtils import *
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


'''
class TuringPath(TuringGraph):
    pass
'''


class TuringGraph:
    def __init__(self,codes,paths=[],**kwargs):
        self.codes=codes
        self.paths=paths
        self.state=cb.BLOCKON
        self.states_register=[]
        self.cid=cb.DEFAULT_CID
        for key in kwargs.keys():
            setattr(self,key,kwargs[key])
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
        


    def add(path):
        #TO DO
        self.paths.append(path)
        pass
    def __str__(self):
        s=""
        for i in self.codes:
            s+=str(i.pos)+":"+str(i.code)+"\n"
        return s
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
        output=bitarray([ True for i in range(2*len(c))])
        
        for i in p:
            c.append(i)
        c.sort()
        j=0
        last_pos=0

        for i in c:
            print i.cid,i
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
                if i.cid==ccid:
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
    paths.append( TuringCode(220,cb.BLOCKOFF))
    paths.append( TuringCode(220,cb.OFF))
    g= TuringGraph(codes)
    p= TuringGraph(paths)

    print g; 
    for i in g.iter_all_path(0):
        print "yield:",i
        l=[]
        for j in i:
            print j
            l.append(TuringCode(g.codes[j].pos,g.codes[j].code))
        print TuringGraph(l)
        print g.translate_path_into_bits(TuringGraph(l))
        a.append(g.translate_path_into_bits(TuringGraph(l)))
        h[g.translate_path_into_bits(TuringGraph(l))]=1


    print g.translate_path_into_bits(p)


    for i,x in enumerate(a):
        for j,y in enumerate(a):
            print i,j
            print x,y
            print twobitarray_and(x,y)
    print h




if __name__=="__main__":
    test()



