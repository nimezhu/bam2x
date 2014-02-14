#!/usr/bin/env python
# Programmer : zhuxp
# Date: 
# Last-modified: 02-13-2014, 19:27:00 EST
import os
import sys
from distutils.core import setup

metadata = {

          'name':"bam2x",
          'version':"0.05",
          'description':"bionformatics python lib for query bam files, xplib final version, mv to bam2x lib in future version",
          'author':"Xiaopeng Zhu",
          'author_email':"nimezhu@gmail.com",
          'packages':[
                    "bam2x",
                    "bam2x.DBI",
                    "bam2x.Run",
                    "bam2x.MRun",
                    "bam2x.Tools",
                    "bam2x.Turing",
                    ],
          'package_dir':{"":"lib"},
          'scripts':[
                   "bin/bam2x"
                   ],
          'install_requires':['numpy>=1.7.0','pysam>=0.7.5','twobitreader>=2.9','bx-python>=0.7.1','bitarray>=0.8.1']


}

def main():
    if not float(sys.version[:3])>2.4:
        sys.stderr.write("CRITICAL: Python version must be greater than or equal to 2.4!")
        sys.exit(1)
    dist=setup(**metadata)




    
if __name__=="__main__":
    main()





