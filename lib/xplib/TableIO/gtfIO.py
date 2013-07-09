#!/usr/bin/python
# Wei Guifeng, <guifengwei@gmail.com>
# Revised By zhuxp <nimezhu@gmail.com>
# -- coding: utf-8 --
'''
Gencode annotation or Cufflinks transcript.gtf

usage:
    for gene in TableIO.parse("cufflinks.transcript.gtf", 'gtf'):
        print "{0}\t{1}\t{2}".format(gene.chr, gene.strand, (gene.stop - gene.start))
        # it will print gene in the format of GenePred(Genetab, e.g GeneBed class)
        print gene
'''
from xplib.Annotation import GeneBed
import SimpleIO
import sys, re, string, types, copy

def GTFIterator(handle,**kwargs):
    ''' gtf line iterator(output format: GenePred) '''
    gene={}
    for line in SimpleIO.SimpleIterator(handle,**kwargs):
            #transid = re.findall(r'transcript_id \"([\w\.]+)\"',line[8])[0].strip()
            if line[2] == 'transcript':
                if gene:
                    yield gtfMakeGeneBed(gene)
                gene={}
                gene['chr'] = line[0]
                gene['start'] = int(line[3]) - 1
                gene['stop']  = int(line[4])
                gene['strand']= line[6]
                gene['id']  = re.findall(r'transcript_id \"([\w\.]+)\"',line[8])[0].strip()
                gene['name2'] = re.findall(r'gene_id \"([\w\.]+)\"',line[8])[0].strip()
                gene['exon']=[]
                gene['cds'] =[]
                continue
            else:
                id_now  = re.findall(r'transcript_id \"([\w\.]+)\"',line[8])[0].strip()
                if not gene.has_key("id") or id_now!=gene['id']:
                    if gene:
                        yield gtfMakeGeneBed(gene)
                    gene={}
                    gene['chr'] = line[0]
                    gene['start'] = int(line[3]) - 1
                    gene['stop']  = int(line[4])
                    gene['strand']= line[6]
                    gene['id']  = re.findall(r'transcript_id \"([\w\.]+)\"',line[8])[0].strip()
                    gene['name2'] = re.findall(r'gene_id \"([\w\.]+)\"',line[8])[0].strip()
                    gene['exon']=[]
                    gene['cds'] =[]
              
            if line[2] == 'exon':
                # gtf: 1-based; bed: 0-based
                gene['exon'].append( [line[0], str(int(line[3])-1), line[4]] )
                if gene['start'] > int(line[3])-1 :
                    gene['start']=int(line[3])-1
                if gene['stop'] < int(line[4]):
                    gene['stop']=int(line[4])
                continue
            if line[2] == 'CDS':
                gene['cds'].append( [str(int(line[3])-1), line[4]] )
    yield gtfMakeGeneBed(gene)

def gtfMakeGeneBed(data):
    ''' '''
    header = ['id', 'chr', 'strand', 'start', 'stop', 'cdsStart', 'cdsEnd', \
              'exonCount','exonStarts', 'exonEnds', 'score', 'name2' ]
    GeneTab_line = copy.copy(data)
    GeneTab_line['score'] = 0
    GeneTab_line['exonCount'] = len( GeneTab_line['exon'] )
    starts, stops =[] , []
    for e in sorted(GeneTab_line['exon'], key=lambda x: int(x[1]) ):
        starts.append( int(e[1]) )
        stops.append( int(e[2]) )
    GeneTab_line['exonStarts'] = ','.join(map(str,starts)) + ','
    GeneTab_line['exonEnds'] = ','.join(map(str, stops)) + ','
    s,t=[],[]
    if GeneTab_line['cds']:
        for i in GeneTab_line['cds']:
            s.append(int(i[0]))
            t.append(int(i[1]))
        GeneTab_line['cdsStart']=min(s)
        GeneTab_line['cdsEnd']  =max(t)
    else:
        GeneTab_line['cdsStart']=stops[-1]
        GeneTab_line['cdsEnd']  =stops[-1]
    
    x=[]
    for g in header:
        x.append(GeneTab_line[g])
    gene = GeneBed(x)
    return gene


