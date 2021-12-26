#! /remote/gpu05/anaconda3/bin/python3
# -*- coding: utf-8 -*-

from bedSep import *
from userFunc import * 

# =============================================================================
#                  define global parameters
# =============================================================================

# path to the .bed file
if len(sys.argv) == 1:
    raise Exception("No input file")
else:
    bedPath = sys.argv[1]

# target chromosome index
chrRng = ['1'] 
# target chromosome length
chrLength = [3188341]

# =============================================================================
#                       testRun begin
# =============================================================================

# separate the .bed file
bedSep(bedPath)

# get radial distribution function and config file
target_file_path = "./"
getGrAndConfig(target_file_path,50000,25000)
os.system("rm 1getGr.out")
os.system("rm 2getConfig.out")

# Run ISS
os.system("python3 ./runISS.py chr1/chr1-0 testRun 1")

# compute compressibility
namesubfix="finres1"
namecompress = "compress1_test.out"

calConformation(namesubfix,150)
calCompressibility(namecompress)
os.system("rm 5getConfor.out")
os.system("rm chr1/*.confor")

# return final file
summarizeAllValue(namesubfix,namecompress,"Potential1_test.like_bed")
      
    
          