import sys,os,argparse
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import spline

def ParseArg():
    p=argparse.ArgumentParser(description = "draw disto/phasogram forprocessed count data",epilog = "Library dependency: matplotlib, scipy, numpy")
    group=p.add_mutually_exclusive_group()
    group.add_argument("-d","--distogram",action='store_true',help="draw distogram for MNase_seq data")
    group.add_argument("-p","--phasogram",action='store_true',help="draw phasogram for MNase-seq data")
    p.add_argument("input",type=str,metavar='input_count',help='input txt file for processed count data')
    p.add_argument("-x","--xlim",type=int, nargs='+',default=[-200,2000],dest="xlim",help="range for x axis: min_x, the left bound; max_x, the right bound. (default: -200,2000)")
    p.add_argument("-o","--output",type=str,dest="output",help="the output figure file format, can be format of emf, eps, pdf, png, ps, raw, rgba, svg, svgz")
    if len(sys.argv)==1:
        print >>sys.stderr,p.print_help()
        exit(0)

    return p.parse_args()

def Main():
    args=ParseArg()
    data=np.loadtxt(args.input,delimiter='\t',dtype=float)
    min_x=int(args.xlim[0])
    max_x=int(args.xlim[1])
    start=int(data[0,0])
    peak=data[:,1].argmax()
    plt.ioff()
    plt.plot(np.array(range(min_x,max_x)),data[(min_x-start):(max_x-start),1],color='r',label="real_count")
    if args.distogram:
        plt.annotate('local max: '+str(peak+start)+"bp",xy=(peak+start,data[peak,1]),xytext=(peak+start+30,0.8*data[peak,1]),)
                #   arrowprops=dict(facecolor='black', shrink=0.05))

    # smoth the plot
    xnew=np.linspace(min_x,max_x,(max_x-min_x)/5)
    smooth=spline(np.array(range(min_x,max_x)),data[(min_x-start):(max_x-start),1],xnew)
    plt.plot(xnew,smooth,color='g',label='smooth(5bp)')
    max_y=max(data[(min_x-start):(max_x-start),1])
    min_y=min(data[(min_x-start):(max_x-start),1])
    plt.xlabel("Distance")
    plt.ylabel("Counts")
    plt.xlim(min_x,max_x)
    plt.ylim(min_y*0.9,max_y*1.1)
    plt.title(os.path.basename(args.input).split("_"+str(start))[0])
    plt.legend()
    plt.savefig(os.path.basename(args.input).split("_"+str(start))[0]+"_%d~%dbp."%(min_x,max_x)+args.output)
    print >>sys.stderr,"output figure file generated!!"

if __name__ == '__main__':
    Main()

