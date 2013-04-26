#!/usr/bin/python
# programmer : zhuxp
# usage: sys.argv[0]  haploblock.tab  gene.tab chr1:1-100000 
'''
 draw figures of haploblock and genes
 input file : gene.tab (download from UCSC table)
              haploblock distribution file ( Generated from ReadRaw.py and *.simple table)
                            *.simple table were generated from xHaploBlockBinomialTest.py (pHaploBlockBinomialTest.sh)
                            the input for that program is BAM file and HapCUT result file
'''

'''
KNOWN ISSUES: No gridspec in sysbio
'''

import sys
import numpy as np
import matplotlib.pyplot as plt

import matplotlib.gridspec as gridspec
from zlab.zbed import *

def addGeneToFig(gene,ax,name=0):
    '''
    add gene to figures
    from green to red is the gene direction
    '''
    if name==1:
        ax.text(gene.start,15,gene.id)
    cds=gene.new_cds()
    utr5=gene.new_utr5()
    utr3=gene.new_utr3()
    cds_exons=cds.Exons()
    for cds_exon in cds_exons:
        ax.bar(cds_exon.start,10,cds_exon.stop-cds_exon.start,color="blue",alpha=0.5)
        ax.bar(cds_exon.start,-10,cds_exon.stop-cds_exon.start,color="blue",alpha=0.5)
    if not utr3 is None:
        for utr3_exon in utr3.Exons():
            ax.bar(utr3_exon.start,5,utr3_exon.stop-utr3_exon.start,color="red",alpha=0.5)
            ax.bar(utr3_exon.start,-5,utr3_exon.stop-utr3_exon.start,color="red",alpha=0.5)
    if not utr5 is None:
        for utr5_exon in utr5.Exons():
            ax.bar(utr5_exon.start,5,utr5_exon.stop-utr5_exon.start,color="green",alpha=0.5)
            ax.bar(utr5_exon.start,-5,utr5_exon.stop-utr5_exon.start,color="green",alpha=0.5)
    for intron in gene.Introns():
        ax.bar(intron.start,1,intron.stop-intron.start,color="yellow",alpha=0.5)
        ax.bar(intron.start,-1,intron.stop-intron.start,color="yellow",alpha=0.5)
        #print intron


def parse(x):
    '''
    parse x:y in Blocks
    '''
    a=x.split(" ")
    b=a[0].split(":")
    c=[0,0,"."]
    c[0]=int(b[0])
    c[1]=int(b[1])
    if len(a)>1:
        c[2]=a[1]
    return c

class tBlock(Bed):
    def __init__(self,x):
        self.chr=x[0].rstrip()
        self.start=int(x[1])
        self.stop=int(x[2])
        self.id=x[3]
        self.copy_numbers=parse(x[4])
        self.hEpi=dict()
        self.hEpi["H3K27ac"]=parse(x[5])
        self.hEpi["H3K36me3"]=parse(x[6])
        self.hEpi["H3K4me1"]=parse(x[7])
        self.hEpi["H3K4me3"]=parse(x[8])
        self.hEpi["H3K9me3"]=parse(x[9])
        self.allele_flag=0
        for i in self.hEpi.keys():
            a=self.hEpi[i]
            if a[2]=="*":
                self.allele_flag=1




def addToFig(block,ax,pos,j,d):
        k=block.hEpi.keys()
        k.sort()
        phase=[list(),list()]
        for i in k:
            a=block.hEpi[i]
            phase[0].append(a[0])
            phase[1].append(-a[1])
        print phase[0]
        print phase[1]
        print k
        ind=np.arange(5)
        ind+=pos
        width=1
        rects0=ax.bar(ind,phase[0],width,color=["r","g","b","y","black"])
        rects1=ax.bar(ind,phase[1],width,color=["r","g","b","y","black"])
        ax.text((j+0.5)/21,0.8,"%d"%d,transform=ax.transAxes) 
        for i,x in enumerate(k):
            a=block.hEpi[x]
            if a[2]=="*" and a[0]>a[1]:
                 height = rects0[i].get_height()
                 ax.text(rects0[i].get_x()+rects0[i].get_width()/2., 1.05*height, '*',ha='center', va='bottom')
            if a[2]=="*" and a[1]>a[0]:
                 height = rects1[i].get_height()
                 ax.text(rects1[i].get_x()+rects0[i].get_width()/2., -1.05*height, '*',ha='center', va='top')
        

        
            
                
        
        
def region_parse(r):
    r=r.strip()
    a=r.split(":")
    b=a[1].split("-")
    return (a[0],int(b[0]),int(b[1]))

def Main():
    region=sys.argv[3]
    (region_chr,region_start,region_stop)=region_parse(sys.argv[3])
    regionBed=Bed([region_chr,region_start-1,region_stop,"QueryRegion",".",0])
    '''
    READING GENE into LIST
    '''
    genes=list()
    genefile=open(sys.argv[2])
    for line in genefile:
        line=line.strip()
        if line[0]=="#": continue
        if len(line)==0: continue
        gene=GeneBed(line.split("\t"))
        if Bed.overlap(gene,regionBed):
            genes.append(gene)
            print gene
    genes.sort()
    '''
    READING BLOCK INTO LIST
    '''
    f=open(sys.argv[1])
    blocks=list()
    for line in f:
        line=line.strip()
        if line[0]=="#": continue
        if len(line)==0: continue
        x=line.split("\t")
        if not x[1].isdigit: continue
        try:
            
            block=tBlock(x)
        except:
            continue
        if Bed.overlap(block,regionBed):
            blocks.append(block)
    blocks.sort()

    
    '''
    initialize fig and draw the first panel of chromosome 
    '''
    number=len(blocks)
    panel_number=number/20+1
    if number%20==0: panel_number-1
    fig = plt.figure(figsize=(25,panel_number*4+3), frameon = True )
    ratios=[4 for i in range(panel_number+2)]
    ratios[0]=2
    ratios[1]=1
    gs = gridspec.GridSpec(panel_number+2, 1, width_ratios=[1],height_ratios=ratios)
    '''
    Draw Blocks Distribution in Chromosome in the first panel
    '''
    ax1=plt.subplot(gs[1])
    ax1.set_yticks([])
    ax1.axis([regionBed.start,regionBed.stop,0,10])
    ax1.set_ylim(0,10)
    ax1.set_xticks([])
    for i,x in enumerate(blocks):
        if x.allele_flag==1:
            ax1.bar(x.start,10,x.stop-x.start,color="r")
        else:
            ax1.bar(x.start,10,x.stop-x.start,color="b")
        if (i+1)%5==0: 
            height=-3
            ax1.text((x.start+x.stop)/2,height,'%d'%(i+1),ha='center',va='bottom')
    '''
    Draw Genes Distribution in the second panel
    '''
    ax2=plt.subplot(gs[0])
    ax2.text(0.05,0.7,region_chr+":"+str(region_start)+"-"+str(region_stop),transform=ax2.transAxes,color="black") 
    ax2.set_yticks([])
    ax2.axis([regionBed.start,regionBed.stop,-30,30])
    for i in genes:
        addGeneToFig(i,ax2,1)
    
    '''
    Draw blocks in other panels 
    '''
    #gsBlocks=gs[2:]
    ax=[0 for i in range(panel_number)]

    for i in range(panel_number):
        ax[i]=plt.subplot(gs[i+2])
    for i,x in enumerate(blocks):
        addToFig(x,ax[i/20],(i%20)*10+5,i%20,i+1)
        ax[i/20].set_xticks([])    
        ax[i/20].set_xlim(0,210)
    
    ax[-1].text(0.95,0.9,"H3K27ac",transform=ax[-1].transAxes,color="r") 
    ax[-1].text(0.95,0.7,"H3K36me3",transform=ax[-1].transAxes,color="g") 
    ax[-1].text(0.95,0.5,"H3K4me1",transform=ax[-1].transAxes,color="b") 
    ax[-1].text(0.95,0.3,"H3K4me3",transform=ax[-1].transAxes,color="y") 
    ax[-1].text(0.95,0.1,"H3K9me3",transform=ax[-1].transAxes,color="black") 
 #   ax[-1].text(0.95,0.9,["H3K27ac","H3K36me3","H3K4me1","H3K4me3","H3K9me3"],transform=ax[-1].transAxes,color=["r","g","b","y","black"]) 
 #   ax[-1].legend()
    '''
    Show Figures
    '''
    plt.show()    
    fig.savefig(region_chr+"_"+str(region_start)+"_"+str(region_stop)+".png")


    
if __name__=="__main__":
    Main()
