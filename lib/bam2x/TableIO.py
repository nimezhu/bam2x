import csv
import bam2x
import sys
from bam2x import IO,Annotation
from bam2x import Translator
hclass = {
    "bed3":Annotation.BED3,
    "bed6":Annotation.BED6,
    "bed12":Annotation.BED12,
    "bed":Annotation.BED6,
    "vcf":Annotation.VCF,
    "bed4":Annotation.BED4,
}
htranslate = {
    "bam2bed12": Translator.BamToBed12,
    "bam2fragment":Translator.BamToFragmentIterator,
}
FormatToIterator=dict(hclass.items()+htranslate.items())
def parse(handle,convert_cls=None,**dict):
    if hclass.has_key(convert_cls):
        return parse_tuples(handle,hclass[convert_cls],**dict)
    elif htranslate.has_key(convert_cls):
        return htranslate[convert_cls](handle,**dict)
    elif convert_cls is not None:
        return parse_tuples(handle,convert_cls,**dict)
    else:
        return parse_simple(handle,**dict)

def parse_simple(handle,**dict):
    sep="\t"
    if dict.has_key("sep"):
        sep=dict["sep"]
    if isinstance(handle,str):
        try:
            handle=IO.fopen(handle,"r")
            for i in csv.reader(handle,delimiter=sep):
                if len(i)==0: continue 
                l=i[0].strip()
                if len(l)>0 and l[0]=="#": continue
                yield tuple(i)
            handle.close()
        except IOError as e:
            print >>sys.stderr,"I/O error({0}): {1}".format(e.errno, e.strerror)
    else:
        try:    
            for i in csv.reader(handle,delimiter=sep):
                if len(i)==0 :continue
                l=i[0].strip()
                if len(l)>0 and l[0]=="#": continue
                yield tuple(i)
        except:
            raise
    
def parse_tuples(handle,cls,**dict):
    if type(cls)==type("string"):
        for i in parse_simple(handle,**dict):
            yield i
    else:    
        for i in parse_simple(handle,**dict):
            yield cls._make(cls._types(i))


def Main():
    for i in parse(sys.argv[1],"bed12"):
        print i
if __name__=="__main__":
    Main()








