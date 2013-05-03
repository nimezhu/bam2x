import numpy as np
from xplib import TableIO
import pysam
from xplib.Annotation import Bed
import sys,os,argparse
import scipy.cluster.hierarchy as sch
import pylab
from scipy.cluster.vq import kmeans,vq
import matplotlib.pyplot as plt
import matplotlib
from collections import Counter

def ParseArg():
    ''' This Function Parse the Argument '''
    p=argparse.ArgumentParser( description = 'plot heatmap or average patterns for intensities around flanking regions of interested locations', epilog='Library dependency : xplib, pylab, scipy, numpy, pysam, matplotlib')
    group=p.add_mutually_exclusive_group()
    group.add_argument("-H","--Heatmap",action='store_true',help="draw heatmap for intensities of all regions")
    group.add_argument("-A","--Average",action='store_true',help="draw average pattern for intensities of all regions")
    p.add_argument('-b','--bams',nargs='+',dest="bams",type=str,help="the list of bam files containing raw reads information (only use the first one if choose -H")
    p.add_argument('-i','--intervals',nargs='+',dest='intervals',type=str,help="the list of bed files containing location information of intervals (only use the first one if choose -H")
    p.add_argument('-d','--direction',action='store_true',help='consider direction of each interval, useful for TSS around regions (default: False)')
    p.add_argument('-r','--resolution',type=int,default=5,help='the resolution for counting reads (default: 5)')
    p.add_argument('-l','--length',type=int,default=1000,help='the length extending from the center to both directions to be drawn (default: 1000)')
    p.add_argument('-f','--frag_l',type=int,default=300,help='the average length of seqeuncing fragments, used to determine center location of each read (default: 300, 150 for MNase-seq)')
    p.add_argument('-n','--name',nargs='+',dest='names',type=str,help="names of the groups (either intervals or bams, can be multiple for only one type) to draw average patterns)")
    p.add_argument('-w','--win_l',type=int,default=3,help='smooth window length for counts in each interval, (default:3, no smooth)')
    p.add_argument('-o','--output',type=str,default='test',help='suffix of output figure file,can be (pdf, eps, png,jpg,...) final will be average_* or heatmap_*')
    if len(sys.argv)==1:
        print >> sys.stderr, p.print_help()
        sys.exit(0)
    return p.parse_args()

def find_center(read,half_len):
    '''
    find center of each fragment from read
    '''
    if not read.is_reverse:
        center=read.pos+half_len
    elif read.is_reverse:
        center=read.aend-half_len
    else:
        print >>sys.stderr,read.id,"No strand"
    return center

def smooth(x,window_len=11,window='hanning'):
        '''
        smooth the count for each interval.
        from:
        http://stackoverflow.com/questions/5515720/python-smooth-time-series-data
        '''
        if x.ndim != 1:
                raise ValueError, "smooth only accepts 1 dimension arrays."
        if x.size < window_len:
                raise ValueError, "Input vector needs to be bigger than window size."
        if window_len<3:
                return x
        if not window in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']:
                raise ValueError, "Window is on of 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'"
        s=np.r_[2*x[0]-x[window_len-1::-1],x,2*x[-1]-x[-1:-window_len:-1]]
        if window == 'flat': #moving average
                w=mp.ones(window_len,'d')
        else:  
                w=eval('np.'+window+'(window_len)')
        y=np.convolve(w/w.sum(),s,mode='same')
        return y[window_len:-window_len+1]


# use center 50bp of the reads
def get_count(interval,bam,r,l,fl,direction=False,win_l=3):
    '''
    generate counts for +-<l> length regions from center of interval
    with <r> resolution. For each read in bam we only use the center
    of corresponding fragment to count (fragment size <fl>)
    '''
    half_len=int(fl/2)
    counts=[]
    for i in interval:
        count=np.zeros(int(2*l/r))
        center=int((i.start+i.stop)/2)
        for j in bam.fetch(i.chr,center-l-(half_len-75),center+l+(half_len+75)):
            j_center=find_center(j,half_len)
            bin=int((j_center-center+l)/r)
            low_bound=min(int(2*l/r),max(0,bin-int(50/r)/2))
            high_bound=max(0,min(int(2*l/r),bin+int(50/r)/2))
            count[low_bound:high_bound]+=1
        if direction and i.strand=='-':
            count=count[::-1]
        elif direction and i.strand!='+':
            print >> sys.stderr, '##ERROR: Strand information not there, please omit -d option'
            sys.exit(0)
        count=smooth(count,window_len=win_l)
        counts.append(count)
        if len(counts)%1000 ==0:
            print >> sys.stderr, "    ##  Counting for interval:",len(counts),"\r",
    print
    return np.array(counts)

global cdict
cdict  = {'red':  ((0.0, 0.0, 0.0),
                   (0.1, 0.0, 0.0),
                   (0.5, 1.0, 1.0),
                   (0.9, 1.0, 1.0),
                   (1.0, 0.8, 0.8)),

         'green': ((0.0, 0.0, 0.0),
                   (0.1, 0.0, 0.0),
                   (0.5, 1.0, 1.0),
                   (0.9, 0.0, 0.0),
                   (1.0, 0.0, 0.0)),

         'blue':  ((0.0, 0.8, 0.8),
                   (0.1, 1.0, 1.0),
                   (0.5, 1.0, 1.0),
                   (0.9, 0.0, 0.0),
                   (1.0, 0.0, 0.0))
        }

def feature_scale(matrix):
    '''
    Feature scaling and standardization
    x=(x-mean(X))/(max(X)-min(X))
    '''
    matrix=(matrix-np.mean(matrix))/(np.max(matrix)-np.min(matrix))
    
    return matrix

def main():
    args=ParseArg()
    # determin if intervals or bams is multiple
    if len(args.bams)>1 and len(args.intervals)>1:
        print "Only one of intervals or bams can have multiple files!!"
        sys.exit(0)
    elif len(args.bams)>1:
        m_signal="bam" # signal for which one is multiple
    elif len(args.intervals)>1:
        m_signal="intervals"
    else:
        m_signal='none'

    #Heatmap option only allow one bam
    if m_signal=='bam' and args.Heatmap:
        print "Heatmap option only allows one bam!!"
        sys.exit(0)

    # store bam files with indexing and count information
    bams={}
    read_numbers={}
    print >> sys.stderr, "##  Starting read and index BAM files: "
    for i in range(len(args.bams)):
        if m_signal=='bam':
            try:
                temp_name=args.names[i]
            except TypeError:
                print >>sys.stderr, "##ERROR: you need to add '-n' for multiple BAM files"
                sys.exit(0)
        else:
            temp_name='bam'
        print >> sys.stderr,"  ##  Indexing for bam file of '"+temp_name+"'"
        bams[temp_name]=pysam.Samfile(args.bams[i],'rb')
        print >> sys.stderr,"    ##  counting total reads number <slow>"
        ss=0
        for chr in bams[temp_name].references:
            ss+=bams[temp_name].count(chr)
        read_numbers[temp_name]=ss
        print >> sys.stderr

    # store interval files
    intervals={}
    print >> sys.stderr, "##  Starting reading intervals:"
    if m_signal=='intervals':
        interval_names=args.names
    else:
        interval_names=['interval']
    if m_signal!='bam':
        for i in range(len(args.intervals)):
            try:
                temp_name=interval_names[i]
            except TypeError:
                print >>sys.stderr, "##ERROR: you need to add '-n' for multiple interval files"
                sys.exit(0)
            print >> sys.stderr,"  ##  Reading for interval file of '"+temp_name+"'\r",
            intervals[temp_name]=TableIO.parse(args.intervals[i],'bed')
            print >> sys.stderr

    resol=args.resolution
    leng=args.length
    frag_l=args.frag_l
    # draw heatmap
    if args.Heatmap:
        print >> sys.stderr,"## Start count reads"
        collect=[]
        interval_n=[0]
        fig=pylab.figure(figsize=(4,8))
        for name in interval_names:
            H_counts=get_count(intervals[name],bams['bam'],resol,leng,frag_l,args.direction,args.win_l)
            H_counts=H_counts*5E7/read_numbers['bam']
            H_counts=np.log(H_counts+1)
            #H_counts=feature_scale(H_counts)
            centroids,_=kmeans(H_counts,5)
            idx,_=vq(H_counts,centroids)
            order=[i[0] for i in sorted(enumerate(idx), key=lambda x:x[1])]
            H_counts=H_counts[order,:]
            collect.append(H_counts)
            interval_n.append(H_counts.shape[0])
        collect=np.vstack(collect)
        print >> sys.stderr,"## Start draw heatmap for intereval"
        #heatmap
        fig=pylab.figure(figsize=(4,8))
        axmatrix=fig.add_axes([0.2,0.1,0.65,0.8])
        plt.register_cmap(name='cust_cmap', data=cdict)
        collect=np.minimum(collect,np.mean(collect)+3*np.std(collect)) # remove outlier with are extremely large
        im=axmatrix.matshow(collect,aspect='auto',origin='lower',cmap=plt.get_cmap('Reds'))
        axmatrix.set_xticks([])
        axmatrix.set_yticks([])
        
        cum_interval_n=np.cumsum(interval_n)
        for n in cum_interval_n:
            axmatrix.axhline(y=n,color='black')
        
        #yticks
        aylabel=fig.add_axes([0.18,0.1,0.0,0.8])
        aylabel.set_yticks(cum_interval_n)
        aylabel.set_xticks([])
        #xticks
        axlabel=fig.add_axes([0.2,0.08,0.65,0.0])
        axlabel.set_xticks([-leng,0,leng])
        axlabel.set_yticks([])
        #plot colorbar
        axcolor=fig.add_axes([0.86,0.1,0.03,0.8])
        pylab.colorbar(im,cax=axcolor)
        fig.savefig('heatmap_'+args.output)

    # draw averge patterns
    if args.Average:

        print >> sys.stderr,"##  Start count reads"
        collect={}
        
        if m_signal=='bam':
            for name in args.names:
                intervals['interval']=TableIO.parse(args.intervals[0],'bed')
                print >> sys.stderr,"  ##  count for bam["+name+"]"
                H_counts=get_count(intervals['interval'],bams[name],resol,leng,frag_l,args.direction)                
                collect[name]=np.sum(H_counts,axis=0)/H_counts.shape[0]*5E7/read_numbers[name]
            names=args.names
        else:
            for name in interval_names:
                print >> sys.stderr,"  ##  count for interval["+name+"]"
                H_counts=get_count(intervals[name],bams['bam'],resol,leng,frag_l,args.direction)
                collect[name]=np.sum(H_counts,axis=0)/H_counts.shape[0]*5E7/read_numbers['bam']
            names=interval_names

        
        fig=plt.figure()
        for i,name in enumerate(names):
            col=matplotlib.cm.Paired((i+1)*1.0/(len(names)+2),1)
            plt.plot(np.array(range(-leng,leng,resol))+resol/2.0,collect[name],color=col)
        pylab.legend(names,loc='upper right')
        plt.xlabel('Distance to center')
        plt.ylabel('Average coverage for 5E7 reads')
        fig.savefig('average_'+args.output)
    


if __name__=="__main__":
    main()









