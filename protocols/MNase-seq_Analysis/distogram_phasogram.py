#-------------------------------------------------------------------------------
# Name:        distogram
# Purpose:
#
# Author:      Pengfei
#
# Created:     31/08/2012
# Copyright:   (c) Pengfei 2012
# Licence:     <your licence>
#-------------------------------------------------------------------------------

import sys,os,argparse
import pysam
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import scipy
from scipy.interpolate import spline

plt.ioff()


#for midpoint plot
def get_center(read,half_len):
    #calculate center
    if read.is_reverse:
        center=read.aend-half_len
    else:
        center=read.pos+half_len
    return center



#for re-sampling of data
def find_peak(x,y,mi,ma):
    sx=x[(x>mi)&(x<ma)]
    sy=y[(x>mi)&(x<ma)]
    return sx[np.argmax(sy)]

def detect_peak(x,y):
    n=np.arange(0,6,dtype="float")
    find_peak(x,y,160,200)
    l=np.arange(0,6,dtype="float")
    for i in range(1,6):
        l[i]=find_peak(x,y,190*i-25,190*i+25)
    n=np.array([n,np.ones(6)])
    slope=np.linalg.lstsq(n.T,l)[0]
    return slope[0]
    
    



def ParseArg():
    p=argparse.ArgumentParser( description = "draw distogram/phasogram for MNase_seq data", epilog="Library dependency: pysam, matplotlib, scipy, numpy")
    group=p.add_mutually_exclusive_group()
    group.add_argument("-d","--distogram",action='store_true',help="draw distogram for MNase_seq data")
    group.add_argument("-p","--phasogram",action='store_true',help="draw phasogram for MNase-seq data")
    group.add_argument("-m","--midpoint",action='store_true',help="draw the distribution of MNase-seq midpoint")
    p.add_argument("input",type=str,metavar='input_rmdup',help='dup removed input bam file for mapped MNase-seq data')
    p.add_argument("-b","--input_o",type=str,dest="input_o",metavar="input_bam",help='original input bam file for mapped MNase-seq data,necessary when pile>1')
    p.add_argument("-o","--output",type=str,dest="output",help="the output figure file, can be format of emf, eps, pdf, png, ps, raw, rgba, svg, svgz")
    p.add_argument("-i","-pile",type=int,dest="pile",default=1,help="conditioning the analysis on sites with #[pile] or more read starts (default=1)")
    #p.add_argument("-l","--lambda",type=float,dest="lambd",default=0.05,help="covariance_factor lambda for KDE (default 0.05)")
    p.add_argument("-x","--xlim",type=int, nargs='+',default=[-200,2000],dest="xlim",help="range for x axis: min_x, the left bound; max_x, the right bound. (default: -200,2000)")
    p.add_argument("-n","--num",type=int,dest="num",default=5000000,help="number of distances included to draw density function (default: 5000000)")
    if len(sys.argv)==1:
        print >>sys.stderr,p.print_help()
        exit(0)

    return p.parse_args()

def Main():
    args=ParseArg()
    min_x=int(args.xlim[0])
    max_x=int(args.xlim[1])
    print  >>sys.stderr,"Reading dup-removed bam file..."
    reads=pysam.Samfile(args.input,'rb')
    if args.pile>1:
        print  >>sys.stderr,"Reading original bam file..."
        reads_o=pysam.Samfile(args.input_o,'rb')
    counts=np.zeros(max_x-min_x)
    num=0
    n=0

    #for resampling purpose
    slopes=[]
    count_temp=np.zeros(max_x-min_x)

    for i in reads:
        n=n+1
        if (not i.is_unmapped):
            pointer=reads.tell()
            chrom=reads.getrname(i.tid)
            start=i.pos+1
            
            # Skip chrX and chrY and chrM            
            if (chrom=='chrX' or chrom=='chrY' or chrom=='chrM'):
                continue
            
            # only forward reads for disto and phasogram
            if (not args.midpoint) and (i.is_reverse):
                continue
            
            if args.midpoint:
                center=get_center(i,74)
                start=center

            # for different number of piles: [args.pile]
            # conditioning the analysis on sites with #[args.pile] or more read starts
            if args.pile>1:
                pile=0
                for k in reads_o.fetch(chrom,start-1,start+1):
                    if (not k.is_reverse) and (k.pos+1==start):
                        pile=pile+1
                if (pile<args.pile) or (pile>50):
                    continue
            #print chrom,start,pile            

            for j in reads.fetch(chrom,start+min_x-75,start+max_x+75):
                if args.distogram and (not j.is_unmapped) and (j.is_reverse):
                    if min_x<=(j.aend-start)<max_x:
                        counts[j.aend-start-min_x]+=1
                        num=num+1
                        if num%100000==0:
                            print >>sys.stderr,"Processing:",num,"distances",chrom,start,"read#",n,"\r",
                elif args.phasogram and (not j.is_unmapped) and (not j.is_reverse):
                    if (min_x<=(j.pos+1-start)<max_x) and (j.qname!=i.qname):
                        counts[j.pos+1-start-min_x]+=1
                        count_temp[j.pos+1-start-min_x]+=1
                        num=num+1
                        if num%100000==0:
                            print >>sys.stderr,"Processing:",num,"distances",chrom,start,"read#",n,"\r",

                        #resampling purpose, find standard deviation of phase distance.
                        if num%500000==0:
                            xnew=np.linspace(min_x,max_x,200)
                            smooth=spline(range(min_x,max_x),count_temp,xnew)
                            slope=detect_peak(xnew,smooth)
                            print >>sys.stderr,"In process of %dth 500000 distances,slope is %f"%(num/500000,slope)
                            count_temp=np.zeros(max_x-min_x)
                            slopes.append(slope)
            
                elif args.midpoint and (not j.is_unmapped):
                    center_j=get_center(j,74)
                    if (min_x<=center_j-center<max_x) and (j.qname!=i.qname):
                        counts[center_j-center-min_x]+=1
                        num=num+1
                        if num%100000==0:
                            print >>sys.stderr,"Processing:",num,"distances",chrom,start,"read#",n,"\r",

            if num>=args.num:
                break
            reads.seek(pointer)
    reads.close()
    if args.pile>1:
        reads_o.close()

    print np.sum(counts)
    print counts[0:10]
    print
    print >>sys.stderr,"start drawing..."
    if args.distogram:
        name="Distogram"
    elif args.phasogram:
        name="Phasogram"
    elif args.midpoint:
        name="Midpoint"

    #plt.hist(distancePool,bins=(max_x-min_x),normed=0,range=(min_x,max_x),alpha=0.4)
    #hist,bins=np.histogram(distancePool,bins=(max_x-min_x),range=(min_x,max_x),density=False)
    np.savetxt(name+"_"+args.output.split(".")[0]+"pile-%d_%d~%dbp.txt"%(args.pile,min_x,max_x),np.column_stack((np.array(range(min_x,max_x),int),counts)),delimiter='\t')
    peak=counts.argmax()
    plt.plot(range(min_x,max_x),counts,color='r')
    if args.distogram:
        plt.annotate('local max: '+str(peak+min_x)+"bp",xy=(peak+min_x,counts[peak]),xytext=(peak+min_x+30,0.8*counts[peak]),)
                #   arrowprops=dict(facecolor='black', shrink=0.05))

    # smoth the plot
    xnew=np.linspace(min_x,max_x,200)
    smooth=spline(range(min_x,max_x),counts,xnew)
    
    #resampling purpose
    if args.phasogram:
        np.savetxt(name+"_"+args.output.split(".")[0]+"phase_distances.txt",slopes,delimiter="\t")
        slope=detect_peak(xnew,smooth)
        print "total slope",slope

    plt.plot(xnew,smooth,color='g')

    plt.xlabel("Length")
    plt.ylabel("Density")
    plt.xlim(min_x,max_x)

    plt.title(name+" of "+os.path.basename(args.input).split(".")[0])
    plt.legend()
    plt.savefig(name+"_pile-%d_"%(args.pile)+args.output)
    print >>sys.stderr,"output figure file generated!!"

if __name__ == '__main__':
    Main()





