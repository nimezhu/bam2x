#!/usr/bin/env python
# Programmer : zhuxp
# Date: 
# Last-modified: 03-31-2014, 11:21:27 EDT
from __future__ import print_function
from distutils.core import setup
from setuptools import setup, find_packages
import codecs
import os
import sys
import re

here = os.path.abspath(os.path.dirname(__file__))

def read(*parts):
    # intentionally *not* adding an encoding option to open
        return codecs.open(os.path.join(here, *parts), 'r').read()

def find_version(*file_paths):
    version_file = read(*file_paths)
    version_match = re.search(r"^__version__ = ['\"]([^'\"]*)['\"]",
                             version_file, re.M)
    if version_match:
        return version_match.group(1)
    raise RuntimeError("Unable to find version string.")


metadata = {

          'name':"bam2x",
          'version':find_version('lib',"bam2x",'__init__.py'),
          'description':"bionformatics python lib for query bam files, this version set matplotlib as optional.",
          'author':"Xiaopeng Zhu",
          'license':'GNU General Public License',
          'url':'http://github.com/nimezhu/bam2x',
          'author_email':"nimezhu@gmail.com",
          'packages':[
                    "bam2x",
                    "bam2x.DBI",
                    "bam2x.Run",
                    "bam2x.Run.Plot",
                    "bam2x.MRun",
                    "bam2x.Tools",
                    "bam2x.Turing",
                    ],
          'package_dir':{"":"lib"},
          'platforms':'any',
          'scripts':[
                   "bin/bam2x",
                   "bin/bam2x_ls",
                   ],
          'package_data':{"":["README.md"]},
          'install_requires':['numpy>=1.6.0','pysam>=0.7.0','twobitreader>=2.9','bx-python>=0.7.0','bitarray>=0.8.0'],
}

def main():
    if not float(sys.version[:3])>2.4:
        sys.stderr.write("CRITICAL: Python version must be greater than or equal to 2.4!")
        sys.exit(1)
    dist=setup(**metadata)




    
if __name__=="__main__":
    main()





