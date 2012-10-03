# Programmer : zhuxp
# Date: 
# Last-modified: 03 Oct 2012 16:49:45
from xplib.Annotation import Bed
from xplib.Annotation import TransUnit
import types
def TransUnitIterator(handle,**kwargs):
    if type(handle)==type("s"):
        try:
            handle=open(handle,"r")
        except:
            raise ValueError("Can't open file %s"%handle)
    TU=TransUnit()
    for line in handle:
        line=line.strip()
        if len(line)==0: continue
        if line[0]=="#": continue
        if line=="//" or line=="// ": 
            yield TU
            TU=TransUnit() #Reset
            continue
        x=line.split("\t")

        if x[1]=="OverlapGene:":
            gene=Bed(x[2:])
            TU.append_overlap_gene(gene)
        elif x[1]=="OverlapFeat:":
            feat=Bed(x[2:])
            TU.append_overlap_feat(feat)
        elif x[1]=="NearbyPromoter:":
            TU.append_promoter(x[2])
        elif x[1]=="Promoter Info:":
            TU.promoterInfo=x[2]
        else:
            TU.processHeader(x)

