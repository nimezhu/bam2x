#!/usr/bin/python
# Programmer : zhuxp
# Date: 
# Last-modified: 20 Sep 2012 13:06:47

import os,sys,argparse
from xplib.Annotation import *
from xplib import TableIO
def ParseArg():
    ''' This Function Parse the Argument '''
    p=argparse.ArgumentParser( description = 'Merge Annotations that have same (chromosome,start,stop,strand) and sort them to print', epilog='Library dependency : xplib')
    p.add_argument('-v','--version',action='version',version='%(prog)s 0.2')
    p.add_argument('-o','--output',dest="output",type=str,default="stdout",help="output file")
    p.add_argument('-a','--annotations',dest="db",action="store",default=[],help="feature annotation files",nargs="+",required=True)
    p.add_argument('-A','--dbformat',dest="db_format",action="store",default=[],help="feature annotation files format",nargs="+",required=True)

    
    if len(sys.argv)==1:
        print >>sys.stderr,p.print_help()
        exit(0)
    return p.parse_args()
def Main():
    global args,out
    args=ParseArg()
    if args.output=="stdout":
        out=sys.stdout
    else:
        try:
            out=open(args.output,"w")
        except IOError:
            print >>sys.stderr,"can't open file ",args.output,"to write. Using stdout instead"
            out=sys.stdout
    db_format=args.db_format
    
    if len(db_format)==1:
        db_format=[db_format[0] for i in range(len(args.db))]
    data={}
    for i,f in enumerate(args.db):
        for item in TableIO.parse(f,db_format[i]):
            flag=0
            for feat in Utils.getOverlapFeatures(item,data):
                if feat==item: #define in Class.__cmp__
                    flag=1
            if not flag:
                Utils.appendIntoBinIndex(data,item)
    data_list=[]
    for i in Utils.iterBinIndex(data):
        data_list.append(i)
    data_list.sort()
    for i in data_list:
        print >>out,i



    
if __name__=="__main__":
    Main()




