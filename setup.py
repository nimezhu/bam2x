#!/usr/bin/env python
# Programmer : zhuxp
# Date: 
# Last-modified: 10-24-2013, 18:27:26 EDT
import os
import sys
from distutils.core import setup
def main():
    if not float(sys.version[:3])>2.4:
        sys.stderr.write("CRITICAL: Python version must be greater than or equal to 2.4!")
        sys.exit(1)
    setup(
          name="xplib",
          version="0.1",
          description="bam2x lib",
          author="Xiaopeng Zhu",
          author_email="nimezhu@gmail.com",
          packages=["xplib",
                    "xplib.Annotation",
                    "xplib.DBI",
                    "xplib.Stats",
                    "xplib.Struct",
                    "xplib.TableIO",
                    "xplib.Tools",
                    ],
          package_dir={"":"lib"},
          scripts=["bin/xQuery.py",
                   "bin/xRead.py",
                   "bin/xCmpGene.py",
                   "bin/xGetSeq.py",
                   "bin/xbams2APS.py",
                   ],
          requires=["pysam","twobitreader","bx"],
          )





    
if __name__=="__main__":
    main()





