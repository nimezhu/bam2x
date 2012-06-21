#!/usr/bin/python
# programmer : zhuxp
# usage:
import sys
from getopt import getopt
import networkx as nx
import matplotlib.pyplot as plt

def show_help():
    print >>sys.stderr,"plotCausalGraph.py -m matrix"
    exit(0)
def Main():
    if len(sys.argv)==1: show_help()

    opts,restlist = getopt(sys.argv[1:],"m:oh",\
                        ["matrix=","other","help"])
    for o, a in opts:
        if o in ("-m","matrix"): M = a
        if o in ("-h","--help"): show_help()
    try:
        f=open(M)
    except:
        print >>sys.stderr,"Can't open file",M
        show_help()
    nodes=[]
    edges=[]
    edges_col=[]
    lines=f.readlines()
    for line in lines:
        line=line.strip()
        if line[0]=="#":continue
        a=line.split("\t")
        nodes.append(a[0])
    j=0
 
    for line in lines:
        line=line.strip()
        if line[0]=="#":continue
        a=line.split("\t")
        for i,x in enumerate(a[1:]):
            if i==j:continue
            x=float(x)
            if x>0.5:
                edges.append((nodes[i],nodes[j]))
                edges_col.append("g")
            elif x<0.5:
                edges.append((nodes[i],nodes[j]))
                edges_col.append("r")
        j+=1


    G=nx.DiGraph()
    G.add_nodes_from(nodes)
    G.add_edges_from(edges)
    nx.draw(G,edge_color=edges_col)
    plt.show()



    
if __name__=="__main__":
    Main()

