#!/usr/bin/env python
# Programmer : zhuxp
# Date: 
# Last-modified: 01-27-2014, 15:15:13 EST
import os
import sys
from distutils.core import setup

metadata = {

          'name':"xplib",
          'version':"0.02",
          'description':"bionformatics python lib for query bam files",
          'author':"Xiaopeng Zhu",
          'author_email':"nimezhu@gmail.com",
          'packages':["xplib",
                    "xplib.Annotation",
                    "xplib.DBI",
                    "xplib.Stats",
                    "xplib.Struct",
                    "xplib.TableIO",
                    "xplib.Tools",
                    "xplib.Turing",
                    "xplib.Tuple",
                    "xplib.TuringTuple"
                    ],
          'package_dir':{"":"lib"},
          'scripts':["bin/xQuery.py",
                   "bin/xRead.py",
                   "bin/xCmpGene.py",
                   "bin/xGetSeq.py",
                   "bin/xbams2APS.py",
                   "bin/pBam20f.py",
                   "scripts/bam2rpkm.py",
                   "scripts/bamInfo.py"
                   ],
          #'requires':['pysam (>=0.7.5)','twobitreader (>=2.9)',],
          'install_requires':['pysam>=0.7.5','twobitreader>=2.9','bx-python>=0.7.1']


}

def main():
    if not float(sys.version[:3])>2.4:
        sys.stderr.write("CRITICAL: Python version must be greater than or equal to 2.4!")
        sys.exit(1)
    dist=setup(**metadata)




    
if __name__=="__main__":
    main()





