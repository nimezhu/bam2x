#!/usr/bin/python
# programmer : zhuxp
# usage:
class HapSimple:
    def __init__(self,x):
        self.chr=x[0]
        self.start=int(x[1])# 0 index
        self.stop=int(x[2])
        self.id=x[3]
        self.coverage=int(x[6])
        self.phase0=int(x[7])
        self.phase1=int(x[8])
        self.alleleA=int(x[9])
        self.alleleB=int(x[10])
        self.p1=float(x[11])
        self.p2=float(x[12])
    def __cmp__(self,other):
        return cmp(self.chr,other.chr) or cmp(self.start,other.start)
def Main():
    hHap=dict()
    lName=["H3K27ac","H3K36me3","H3K4me1","H3K4me3","H3K9me3"]
    fName=[]
    print "chr\tstart\tstop\tid\tcopynumbers\t","\t".join(lName)
    for i in lName:
        fName.append(i+".HapBlockBinomialTest.raw.tab")
    for i,fname in enumerate(fName):
        f=open(fname)
        for line in f:
            line=line.strip()
            if line[0]=="#": continue
            if len(line)==0: continue
            x=line.split("\t")
            a=HapSimple(x)
            if hHap.has_key(a.id):
                hHap[a.id][lName[i]]=a
            else:
                hHap[a.id]=dict()
                hHap[a.id][lName[i]]=a
    ids=hHap.keys()
    ids.sort(key=lambda item:hHap[item]["H3K27ac"] )

    for id in ids:
        a=hHap[id]
        print a["H3K27ac"].chr,"\t",a["H3K27ac"].start,"\t",a["H3K27ac"].stop,"\t",id,"\t",
        print str(a["H3K27ac"].alleleA)+":"+str(a["H3K27ac"].alleleB),
        flag=0
        for i in lName:
            print "\t",str(a[i].phase0)+":"+str(a[i].phase1),
            if (a[i].p1<0.1 and a[i].p2<0.1):
                if flag==0:
                    print "*",
                    flag=cmp(a[i].phase0,a[i].phase1)
                else:
                    if flag==cmp(a[i].phase0,a[i].phase1):
                        print "*",
                    else:
                        print "-*",

                #if(a[i].phase0 >= a[i].phase1): print "+",
                #else : print "-",
        print


if __name__=="__main__":
    Main()
