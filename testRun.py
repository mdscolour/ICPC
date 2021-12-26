"""
This is testRun.py file. 
To test the program with the example data, please run:
$ python3 testRun.py exampleData/GSM1542419_chr1.bed

Normally this will takes about 90 minutes. 

It will generate a potential*.like_bed file as final result.
Inside of this file, there are two lines of header and
the third line looks similar with: 
chr1	12500	75000	0	151.000000	3.000000	11.000000	8.000000	0.066213
"""

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
bedSep(chrRng,bedPath)

# get radial distribution function and config file
target_file_path = "./"
getGrAndConfig(chrRng,chrLength,target_file_path,50000,25000)
os.system("rm 1getGr.out")
os.system("rm 2getConfig.out")

# Run ISS for section 0 for chr. 1
os.system("python3 ./runISS.py chr1/chr1-0 testRun 1")

# compute compressibility
namesubfix="finres1"
namecompress = "compress1_test.out"

calConformation(chrRng,namesubfix,150)
print("Creation of conformations success.")

calCompressibility(chrRng,namecompress)
print("Calculation of compressibilities success.")

os.system("rm 5getConfor.out")
os.system("rm chr1/*.confor")

# return final file
namepotres = "Potential1_test.like_bed"
summarizeAllValue(chrRng,namesubfix,namecompress,namepotres)
print("Program finished. Final output name: %s"%namepotres) 
    
          