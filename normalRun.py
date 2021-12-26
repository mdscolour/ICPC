"""
This is normalRun.py file. 
To run the program with the example data, please run:
$ python3 normalRun.py exampleData/GSM1542419_chr1.bed

Directly run this code will takes about 120 hours.

And this is only for the first section in chr. 1.

Change the blocked code in the part "Run ISS" according 
to specific device or cluster to allow parallel computing 
for multiply sections to make the calculation of the 
whole genome affordable.

Further time reduction is possible by changing the noisy 
optimization part in the file "runISS.py". 

At the end it will generate a potential*.like_bed file as 
final result.
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
chrRng = ['1','2','3','4','5','6','7','R'] 
# target chromosome length
chrLength = [3188341,2231883,1799298,1603259,1190845,1033292,949511,2286237]

# =============================================================================
#                       normalRun begin
# =============================================================================

# separate the .bed file
bedSep(chrRng,bedPath)

# get radial distribution function and config file
target_file_path = "./"
getGrAndConfig(chrRng,chrLength,target_file_path,50000,25000)
os.system("rm 1getGr.out")
os.system("rm 2getConfig.out")

# Run ISS
# =============================================================================
os.system("python3 ./runISS.py chr1/chr1-0 normalRun 1")
os.system("python3 ./runISS.py chr1/chr1-0 normalRun 0.1")
# =============================================================================

# compute compressibility
namesubfix="finres0.1"
namecompress = "compress01_100000.out"

calConformation(namesubfix,100000)
print("Creation of conformations success.")

calCompressibility(chrRng,namecompress)
print("Calculation of compressibilities success.")

os.system("rm 5getConfor.out")
os.system("rm */*.confor")

# return final file
namepotres = "Potential01_100000.like_bed"
summarizeAllValue(chrRng,namesubfix,namecompress,namepotres)
print("Program finished. Final output name: %s"%namepotres) 
          