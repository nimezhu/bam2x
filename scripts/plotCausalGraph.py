#!/usr/bin/python 
# programmer : zhuxp
# usage:
import sys
from getopt import getopt
import networkx as nx
import matplotlib.pyplot as plt
def show_help():
    print >>sys.stderr,"\n\nplotCausalGraph.py:  drawing causal graph from LiNGAM output equation matrix B"
    print >>sys.stderr,"Library Dependence:  networkx , matplotlib\n\n"
    print >>sys.stderr,"Usage: plotCausalGraph.py -m matrix -t threshold[default:0.5]"
    print >>sys.stderr,"Options:"
    print >>sys.stderr,"   --matrix,-m file.mat        B matrix file generated from LiNGAM"
    print >>sys.stderr,"   --threshold,-t threshold    draw the line if the absolute value of line weight is bigger than threshold"
    print >>sys.stderr,"   --direct,-d                 draw the arrow"
    print >>sys.stderr,"Matrix File Example:"
    print >>sys.stderr,"H3K9me3\t0"
    print >>sys.stderr,"H3K4me3\t-1\t0"
    print >>sys.stderr,"H3K4me1\t-0.3\t1\t0"
    exit(0)
def Main():
    if len(sys.argv)==1: show_help()

    opts,restlist = getopt(sys.argv[1:],"m:oht:d",\
                        ["matrix=","threshold=","help","direct"])
    threshold=0.5
    direct=False
    for o, a in opts:
        if o in ("-m","--matrix"): M = a
        if o in ("-h","--help"): show_help()
        if o in ("-t","--threshold"): threshold=float(a)
        if o in ("-d","--direct"): direct=True
    if not 'M' in dir():
        show_help()
    try:
        f=open(M)
    except:
        print >>sys.stderr,"Can't open file",M
        show_help()

    max=0
    nodes=[]
    edges=[]
    pos={}
    edges_col=[]
    col={}
    
    rank={}
    lines=f.readlines()
    i=0
    maxcols=0
    for line in lines:
        line=line.strip()
        if line[0]=="#":continue
        a=line.split("\t")
        nodes.append(a[0])
        rank[a[0]]=0
        for k,x in enumerate(a[1:]):
            if k==i: continue
            x=float(x)
            if x>threshold or x<-threshold:
                if k < len(nodes) and rank[a[0]] < rank[nodes[k]]+1:
                    rank[a[0]]=rank[nodes[k]]+1
        if col.has_key(rank[a[0]]):
            col[rank[a[0]]]+=1
        else:
            col[rank[a[0]]]=1
        #pos[a[0]]=(rank[a[0]],col[rank[a[0]]]+rank[a[0]]%2*0.5+rank[a[0]]*0.111)
        pos[a[0]]=(rank[a[0]],col[rank[a[0]]])
        i+=1
    for e in col.values():
        if e > maxcols: maxcols=e
    for e in pos.keys():
        a=pos[e]
        pos[e]=(a[0],float(a[1]+0.05*rank[e]%3)/float(col[rank[e]]+1)*maxcols)





    j=0
 
    G=nx.DiGraph()
    G.add_nodes_from(nodes)
    edges=[]
    for line in lines:
        line=line.strip()
        if line[0]=="#":continue
        a=line.split("\t")
        for i,x in enumerate(a[1:]):
            if i==j:continue
            x=float(x)
            if max<abs(x):max=abs(x)
            if x>threshold:
                edges.append((nodes[i],nodes[j],{'color':'green','weight':x}))
        #        edges_col.append(x)
            elif x<-threshold:
                edges.append((nodes[i],nodes[j],{'color':'red','weight':x}))
        #        edges_col.append(x)
        j+=1
    G.add_edges_from(edges)
    e=G.edges()
    for i in e:
        edges_col.append(G[i[0]][i[1]]['weight'])

    nx.draw(G,edge_cmap=plt.get_cmap("RdYlGn"),edgelist=e,edge_color=edges_col,pos=pos,node_color="y",edge_vmin=-max,edge_vmax=max,linewidths=0,width=2,arrows=direct,node_size=100,font_size=10)
    #nx.draw(G,edge_cmap=plt.get_cmap("RdYlGn"),edge_list=edges,edge_color=edges_col,pos=pos,node_color="y",edge_vmin=-max,edge_vmax=max,linewidths=0,width=2)
    #nx.draw(G,pos=pos,node_color="y",linewidths=0,width=2)
    plt.colorbar()
    plt.show()




    
if __name__=="__main__":
    Main()

