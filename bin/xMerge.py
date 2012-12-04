#!/usr/bin/env python
# Programmer : zhuxp
# Date: 
# Last-modified: 03 Dec 2012 21:36:49
VERSION="0.3"
import os,sys,argparse
from xplib.Annotation import *
from xplib import TableIO
from xplib.Struct import binindex
from xplib import DBI
import time
def ParseArg():
    ''' This Function Parse the Argument '''
    p=argparse.ArgumentParser( description = 'Merge Annotations that have same (chromosome,start,stop,strand) and sort them to print', epilog='Library dependency : xplib')
    p.add_argument('-v','--version',action='version',version='%(prog)s '+VERSION)
    p.add_argument('-o','--output',dest="output",type=str,default="stdout",help="output file")
    p.add_argument('-a','--annotations',dest="db",action="store",default=[],help="feature annotation files",nargs="+",required=True)
    p.add_argument('-A','--dbformat',dest="db_format",action="store",default=[],help="feature annotation files format",nargs="+",required=True)
    p.add_argument('-p','--position',dest="pos",action="store_true",default=False,help="only cosidering the postion, if two annotation are in same postion, treat them as equal")

    
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
    argv=sys.argv
    argv[0]=argv[0].split("/")[-1]
    print >>out,"# This data was generated by program ",argv[0],"(version %s)"%VERSION,
    print >>out,"in bam2x ( https://github.com/nimezhu/bam2x )"
    print >>out,"# Date: ",time.asctime()
    print >>out,"# The command line is :\n#\t"," ".join(argv)
    db_format=args.db_format
    if len(db_format)==1:
        db_format=[db_format[0] for i in range(len(args.db))]
    data=binindex()
    for i,f in enumerate(args.db):
        for item in TableIO.parse(f,db_format[i]):
            flag=0
            for feat in data.query(item):
                if args.pos:
                    if feat.start==item.start and feat.stop==item.stop:
                        flag=1
                elif feat==item: #define in Class.__cmp__
                    flag=1
            if not flag:
                data.append(item)
    data_list=[]
    for i in data:
        data_list.append(i)
    data_list.sort()
    for i in data_list:
        print >>out,i

    
if __name__=="__main__":
    Main()

