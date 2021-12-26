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
chrRng = ['1','2','3','4','5','6','7','R'] 
# target chromosome length
chrLength = [3188341,2231883,1799298,1603259,1190845,1033292,949511,2286237]

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
os.system("python3 ./runISS.py chr1/chr1-0 normalRun 1")
os.system("python3 ./runISS.py chr1/chr1-0 normalRun 0.1")

# compute compressibility
namesubfix="finres0.1"
namecompress = "compress01_100000.out"

calConformation(namesubfix,100000)
calCompressibility(namecompress)
os.system("rm 5getConfor.out")

# return final file
summarizeAllValue(namesubfix,namecompress,"Potential01_100000.like_bed")
      
    
          