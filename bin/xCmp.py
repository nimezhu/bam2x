#!/usr/bin/python
# Programmer : zhuxp
# Date: 
# Last-modified: 13 Sep 2012 20:29:32

import os,sys,argparse
from xplib.Annotation import *
from xplib import TableIO
def ParseArg():
    ''' This Function Parse the Argument '''
    p=argparse.ArgumentParser( description = 'Example: %(prog)s -i [query.bed|query.genebed|query.OddsRatio.out] --input-format [bed|genebed|oddsratiosnp] -a file1.genebed file2.genebed --db_format genebed', epilog='Library dependency : xplib')
    p.add_argument('-v','--version',action='version',version='%(prog)s 0.2')
    p.add_argument('-i','--input',dest="input",type=str,help="input file")
    p.add_argument('-I','--input_format',dest="input_format",action="store",default="bed",help="input file format")
    p.add_argument('-o','--output',dest="output",type=str,default="stdout",help="output file")
    p.add_argument('-a','--annotations',dest="db",action="store",default=[],help="feature annotation files",nargs="+")
    p.add_argument('-A','--annotations_format',dest="db_format",action="store",default=[],help="annotation files format",nargs="+")
    p.add_argument('-m',dest="m",action="store_true",default=False,help="print in table format")

    
    if len(sys.argv)==1:
        print >>sys.stderr,p.print_help()
        exit(0)
    return p.parse_args()
def Main():
    global args,out
    h={} #count the hits number
    hcode={} #count the hits number
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
    if len(db_format)==0:
        db_format=["bed" for i in range(len(args.db))]
    if len(db_format)==1:
        db_format=[db_format[0] for i in range(len(args.db))]
    if len(db_format)!=len(args.db):
        print >>sys.stderr,"the number of annotation files is not same with the number of annotation formats"
        print >>sys.stderr,"db format",db_format
        print >>sys.stderr,"db ",args.db
        exit(0)
    datas=[]
    input_handle=TableIO.parse(args.input,args.input_format)
    beds=[]
    h[args.input]=0
    for i in input_handle:
        beds.append(i)
        h[args.input]+=1
       
    for i,f in enumerate(args.db):
        print >>sys.stderr,"reading",f
        h[f]=0
        fi=TableIO.parse(f,db_format[i])
        '''
        data tuple structure
        (filename,data_binindex)
        usage example:
            bed=Bed([chr,start,end])
            for (file,data) in data:
                for feature in Utils.iterOverlapFeature(bed,data):
                    print feature
        '''
        datas.append((f,Utils.readIntoBinIndex(fi)))
    print >>out,"# Input:",args.input
    for i,f in enumerate(args.db):
        print >>out,"# Feature No."+str(i+1),":",f

    for i in beds:

        print >>out,i,
        code=""
        if not args.m:
            print >>out,""
            print >>out,""
        for (f,d) in datas:
            s=""
            mark=0
            for feature in Utils.iterOverlapFeature(i,d):
                mark+=1
                s+=str(feature)+"\n"
            if mark>0:
                h[f]+=1
                code+="1"
                if not args.m:
                    print >>out,"OVERLAP FEATURES IN",f
                    print >>out,s
                else:
                    print >>out,"\t",1,
            else:
                code+="0"
                if not args.m:
                    print >>out,"NO OVERLAP FEATURES IN",f
                else:
                    print >>out,"\t",0,
        if hcode.has_key(code):
            hcode[code]+=1
        else:
            hcode[code]=1
        if args.m:
            print >>out,"\t","#"+code
        else:
            print >>out,"//"
    for i,f in enumerate(args.db):
        print >>out,"# ",h[f],"/",h[args.input],"overlap with No."+str(i+1),f
    for code in sorted(hcode.iterkeys()):
        print >>out,"# code:",code,"\t",hcode[code]

    
if __name__=="__main__":
    Main()




