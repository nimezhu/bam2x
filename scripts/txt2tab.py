#!/usr/bin/python
# programmer : zhuxp
# usage:
import sys
from getopt import getopt
def show_help():
    print >>sys.stderr,"Usage: txt2tab.py -g genome --bin binsize --txt txtListFile -n "
    exit(0)
def parseGenomeSizes(fn):
    try:
        f=open(fn)
    except:
        print >>sys.stderr,"Can't open genome sizes file",fn
        show_help()
    chr=[]
    length=[]
    i=0
    for line in f:
        line=line.strip()
        a=line.split("\t")
        chr.append(a[0])
        length.append(long(a[1]))
        hChrToCid[a[0]]=i
        i+=1
    genomeSizes['chr']=chr
    genomeSizes['length']=length
    
def init_table():
    '''
    Table[txt_id][chr_id][bin_id]
    '''
    global Table
    Table=[ [] for l in range(len(txtList))]
    for i,chr in enumerate(genomeSizes['chr']):
        bins=genomeSizes['length'][i]/BinSize
        if genomeSizes['length'][i]%BinSize!=0: bins+=1
        for j in range(len(txtList)):
            Table[j].append( [ 0.0 for row in range(bins)])
def read_into_table(fn,txt_id):
    try:
        f=open(fn)
    except:
        print >>sys.stderr,"can't open file", fn
        return 
        #exit(1)
    i=0
    print >>sys.stderr,"Reading No."+str(txt_id)+" File:",fn
    for line in f:
        line=line.strip()
        x=line.split("\t")
        cid=hChrToCid[x[2]]
        pos=long(x[3])
        bin_id=pos/BinSize

        Table[txt_id][cid][bin_id]+=1
        i+=1
        if i%100000==0:
            print >>sys.stderr,"reading",i,"reads\r",
    FileTotalReads[txtList[txt_id]]=i

def normalize_table():
    '''
    Calculate RPKM
    '''
    for txt_id in range(len(Table)):
        total=FileTotalReads[txtList[txt_id]]
        for cid in range(len(Table[txt_id])):
            for bin_id in range(len(Table[txt_id][cid])):
                Table[txt_id][cid][bin_id]*=10000000.0*1000/BinSize/total
                     
    
def print_table():
    print "# Normalize:",Normalize
    for bid in range(len(txtList)):
        print "# File No.",bid,":",txtList[bid],"\t",FileTotalReads[txtList[bid]]
    for cid in range(len(Table[0])):
        for bin_id in range(len(Table[0][cid])):
            print genomeSizes['chr'][cid],"\t",bin_id*BinSize,
            for txt_id in range(len(Table)):
                if Normalize:
                    print "\t%.2f"%Table[txt_id][cid][bin_id],
                else:
                    print "\t%.0f"%Table[txt_id][cid][bin_id],
            print 

def Main():
    global Table,txtList,genomeSizes,FileTotalReads,BinSize,hChrToCid, Normalize
    hChrToCid={}
    genomeSizes={}
    FileTotalReads={}
    BinSize=200
    txtList=[]
    Normalize=False
    if len(sys.argv)==1: show_help()
    opts,restlist = getopt(sys.argv[1:],"b:ohg:n",\
                        ["txt=","genome=","help","bin=","normalize"])
    for o, a in opts:
        if o in ("-b","--txt"): 
            f=open(a)
            for line in f:
                line=line.strip()
                t=line.split("\t")
                txtList.append(t[0])
        if o in ("-g","--genome"): 
            parseGenomeSizes(a)
        if o in ("-h","--help"): show_help()
        if o in ("--bin"): BinSize=int(a)
        if o in ("-n","--normalize"): Normalize=True

    for i in restlist:
        a=i.split('.')
        if a[-1]=="txt":
            txtList.append(i)
    init_table()
    for bid,f in enumerate(txtList):
        read_into_table(f,bid)
    if Normalize:
        normalize_table()

    print_table()




    
if __name__=="__main__":
    Main()

