#!/usr/bin/env pythON
# Programmer : zhuxp
# Date: 
# Last-modified: 11-13-2013, 13:41:04 EST
import xplib.Turing.TuringCodeBook as cb
from bitarray import bitarray 
from xplib.Turing.TuringUtils import *
from xplib.Turing import TuringCode,TuringGraph
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



