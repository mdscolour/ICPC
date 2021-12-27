# Intra-Chromosomal Potential Calculator (ICPC)
This is the program to extract the intra-chromosomal potentials from nucleosomal positioning data.
The obtained results are effective potentials for single nucleosomes, which represent the organization mechanism for that specific section and can serve as a status identifier for further study.
## Related publication
[https://arxiv.org/abs/2112.11785](http://www.picb.ac.cn/hanlab/iNPS.html)
## Dependencies
NumPy, SciPy

Some optional function may need pandas
## Usage
### Test Run
A test run with example data:
```
python3 testRun.py exampleData/GSM1542419_chr1.bed
```
Normally this will take about 90 minutes. 

A success run will generate a potential*.like_bed file as final result.
### Normal Run
A normal run:
```
python3 normalRun.py bed_file_name.bed
```
Directly run this will take about 120 hours.

And this is only for the first section in chr. 1.

Change the blocked code in the part "Run ISS" in file "normalRun.py" 
according to specific device or cluster to allow parallel computing 
for multiply sections to make the calculation of the whole 
genome affordable.

Further time reduction is possible by changing the "noisy 
optimization part" in the file "runISS.py". 

At the end it will generate a potential*.like_bed file as 
final result.
## Algorithm Structure:
The algorithm is organized in the following steps:
* Preparing [iNPS](https://www.nature.com/articles/ncomms5909) data ([iNPS program is here](http://www.picb.ac.cn/hanlab/iNPS.html))
* Obtaining RDF data
* Running simulations to obtain potential parameters
* Computing configurations and compressibilities

In the paper there is a further step:
* Classification

While code of the classification is not included here
but it is already well explained in the paper.
