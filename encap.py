#! /remote/gpu05/anaconda3/bin/python3
# -*- coding: utf-8 -*-

import os
import sys
import copy
import time
import numpy as np
from scipy.interpolate import interp1d
from matplotlib import pyplot as plt
from matplotlib.collections import BrokenBarHCollection
import scipy.cluster.hierarchy as shc
from sklearn.cluster import AgglomerativeClustering
import matplotlib as mpl 
#from pylab import *

#######preparation of .like_bed file, .bed file is needed
def quickBedSep():
    os.system("./0bedSep.py") 
    
def getGrAndConfig(chrRngT,target_file_pathT,seclen=50000,inclen=25000):
    for Chromosome in chrRngT:   
        ### set maximum length
        if Chromosome == '1':	
            maxlen=3188341
        if Chromosome == '2':
            maxlen=2231883
        if Chromosome == '3':	
            maxlen=1799298
        if Chromosome == '4':	
            maxlen=1603259
        if Chromosome == '5':	
            maxlen=1190845
        if Chromosome == '6':	
            maxlen=1033292
        if Chromosome == '7':	
            maxlen=949511
        if Chromosome == 'R':	
            maxlen=2286237
        
        maxsection = int(maxlen/seclen);
        if maxsection*seclen+inclen > maxlen:
            maxsection-=1;
       
       
        ### preparation of .lowconfig and .midgr file

        #if True:        
        ### create .midgrpre and .preconfig. Run in the main folder.
        target_file = "%s/_chr%s.like_bed"%(target_file_pathT,Chromosome)
        ### compile
        os.system("g++ 1getGr.cpp")
        ### create folders if not exist
        try:
            os.system("mkdir chr%s"%Chromosome)
        except:
            pass
        ### run with 6 paramters
        ### the default Gr is calculated to 1000 with stepsize 5
        ### i.e. in range(0,1000,5). Stepsize 5 is according to data quality.
        ### 1000 is decided by the scale of interest. The change of 
        ### these two paremeters need to go into 1getGr.cpp in function getCanGr() 
        os.system("./a.out %s %d %s chr%s %d %d"%
        (Chromosome,maxsection,target_file,Chromosome,seclen,inclen))
        
        #if False:
        ### from .preconfig to .lowconfig. Run in the subfolder.
        os.system("g++ 2getConfig.cpp")
        os.chdir("chr%s"%Chromosome)
        for ita in range(maxsection):
            os.system("../a.out %d %s"%(ita,Chromosome))
        #os.chdir("..")
        
        #if False:
        ### here is step 3
        ### from .midgrpre to .midgr. Run in the subfolder.
        #os.chdir("chr%s"%Chromosome)
        patrng = np.asarray([f for f in os.listdir() if f.endswith("midgrpre")])
        patrng.sort()
    
        res = []
        ### final Gr range decides here
        xres = np.arange(0,600,0.01)## first peak <200, so choose maximum length 600 here
        for pat in patrng:
            x,y = np.genfromtxt(pat).T
            resinp = interp1d(x, y, kind='linear')
            #resinp = interp1d(x, y, kind='quadratic')
            yres = resinp(xres)
            np.savetxt(pat[:-3], np.vstack((xres,yres)).T,fmt='%f')##save deleting "pre"
        
        os.system("rm *pre*")
        os.chdir("..")

### submitting tasks to the ITP clusters
### to run and get the potential, the .finres file.
### changing computational accuracy need to change "candidates" in file "runISS.py".
def submitquest(Chromosome,secrng): 
    os.chdir("chr%s"%Chromosome)
    for i in secrng:
        #dat = np.genfromtxt("/scratch/li/MolecularMC/v8ES1/res/chr2-%i.finres"%i)
        #x0 = dat[-1,1:-1]
        #snx0 = "%i %i %i %i"%(x0[0],x0[1],x0[2],x0[3])
        
        os.system("cp ../opt_ncxxx_itp.sh ./ncchr%s-%d.sh"%(Chromosome,i))
        os.system('sed -i "s/xxx/chr%s-%d/g" ./ncchr%s-%d.sh'%(Chromosome,i,Chromosome,i))
        os.system('sed -i "s/yyy/chr%s/g" ./ncchr%s-%d.sh'%(Chromosome,Chromosome,i))
        #exit()
        os.system("qsub ./ncchr%s-%d.sh"%(Chromosome,i))
        os.system("rm ./ncchr%s-%d.sh"%(Chromosome,i)) 
    os.chdir("..")

### submit all section that is not computed. 
### search according to .midgr files and .namesubfix files
def submitAllMissing(namesubfix="finres"):
    chrRng = ['1','2','3','4','5','6','7','R'] 
    for ichr in chrRng:
        os.chdir("chr%s"%ichr)
        isecrng = [int(f[5:-6]) for f in os.listdir() if f.endswith(".midgr")]
        #print(np.amax(isecrng))
        misarr = []
        for i in range(np.amax(isecrng)+1):
            try:
                data = np.genfromtxt("chr%s-%d.%s"%(ichr,i,namesubfix))
                #print(i,data)
            except:
                #print(i," is missing.")
                misarr.append(i)
        print(ichr,misarr)
        os.chdir("..")
        
        submitquest(ichr,misarr)

### calculate conformations for all sections
### it will temporarily change the result file name
def calConformation(namesubfix="finres"):
    os.system("g++ 5getConfor.cpp")
    #chrRng = ['R'] 
    #chrRng = ['1','2','3','4','5','6','7','R'] 
    global chrRng
    
    for ichr in chrRng:
        os.chdir("chr%s"%ichr)
        isecrng = [int(f[5:-6]) for f in os.listdir() if f.endswith(".midgr")]
        #isecrng = [int(f[5:-9]) for f in os.listdir() if f.endswith(".finres01")]
        isecrng.sort()
        #print(isecrng)
        for i in isecrng:
            ### third parameter is the conformation number, times 10 is MSC
            #os.system("mv chr%s-%d.%s chr%s-%d.finrestmp"%(ichr,i,namesubfix,ichr,i))
            os.system("../a.out %d %s 100000 confor %s"%(i,ichr,namesubfix))  
            #os.system("mv chr%s-%d.finrestmp chr%s-%d.%s"%(ichr,i,ichr,i,namesubfix)) 
        os.chdir("..")

def calCompress(ichr,ita,boxrngT):
    dataname = "chr%s-%d.confor"%(ichr,ita)
    
    secSizeStart = np.nan
    secSizeEnd = np.nan
    with open("chr%s-%d.lowconfig"%(ichr,ita), "r") as infile:
        for iconfig in infile:
            data = iconfig.split(":")
            if data[0]=="global_boxl":
                secSizeStart = float(data[1])
            if data[0]=="global_boxr":
                secSizeEnd = float(data[1])
                break
    sectionsize = int(secSizeEnd-secSizeStart)
    if np.isnan(sectionsize):
        print("error reading config file in calCompress.")
        exit()
        
    #stdarr = []
    rhoktKTarr = []
    for boxsize in boxrngT:
        #binshist = np.arange(0,75000-1+boxsize,boxsize)
        binshist = np.linspace(0,sectionsize,boxsize+1)
        #binshistmid = (binshist[:-1]+binshist[1:])/2.0
        density = [] 
        with open(dataname)as f:            
            for line in f:
                datline = np.asarray(line.split()).astype(float)
                nsam = len(datline)
                ### counting
                ht = np.histogram(datline,binshist)
                density.append(ht[0])
        #print(boxsize,np.array(density).shape,np.sum(density[0]),density[0][0])
        denbinshist = np.arange(np.amin(density)-0.5,np.amax(density)+0.6,1
            )*boxsize*1.0/nsam
        ### real density array
        density = np.asarray(density).flatten()#*boxsize*1.0/75000
        #ktKT = 75000*np.var(density,ddof=1)/np.mean(density)#**2
        rhoktKT = np.var(density)/np.mean(density)
        rhoktKTarr.append(rhoktKT)
        #stdarr.append(np.std(density,ddof=1))
        #print(boxsize,len(density),stdarr[-1])
#            ht = np.histogram(density*75000/nsam,denbinshist,density=True)
#            plot((denbinshist[:-1]+denbinshist[1:])/2.0,ht[0],".-",label=r"$M_b$=%d"%boxsize)
#        xlabel(r"$\rho/\rho_0$")
#        ylabel(r"distribution function $P(\rho/\rho_0)$")
#        legend()
    #figure()
    #plot(75000/boxrngT,stdarr,'*-',color=col[i])
    x = boxrngT#/75000.
    y = np.array(rhoktKTarr)
    y = y[x>15]
    x = x[x>15]        
    A = np.vstack([x, np.ones(len(x))]).T
    sol = np.linalg.lstsq(A, y,rcond=None)
    v, c = sol[0]
    #print(i,v,c)
    return [ita,v,c]+rhoktKTarr
def calCompressibility(nameout="compress.out"):
    #boxrng = np.arange(1,41,1)
    boxrng = np.arange(16,41,2)
    
    global chrRng
    #chrRng = ['1','2','3','4','5','6','7','R'] 
    #chrRng = ['1'] 
    for ichr in chrRng:
        os.chdir("chr%s"%ichr)
        isecrng = [int(f[5:-6]) for f in os.listdir() if f.endswith(".midgr")]
        #isecrng = [int(f[5:-9]) for f in os.listdir() if f.endswith(".finres01")]
        isecrng.sort()
        #print(isecrng)
        
        ffin = open(nameout,"w")
        ffin.write("0 0 0 ")
        for item in boxrng:
            ffin.write("%d "%item)
        ffin.write("\n")
            
        for i in isecrng:
            #print(ichr,i)
            irescomp = calCompress(ichr,i,boxrng)
            print(ichr,irescomp[:3])
            for itemires in irescomp:
                ffin.write("%f "%itemires)
            ffin.write("\n")
        ffin.close() 
        os.chdir("..")
        
    
if __name__ == '__main__':
    ### get iNPS data, here is only a function for chromosome separation
    #quickBedSep() ### separate _chrX.like_bed from .bed file
    
    
    
    ### get radial distribution function and config file   
    ### last two parameter decide the section length
    ### change section length need to change here
    #chrRng = ['1','2','3','4','5','6','7','R']
    #target_file_path = "/scratch/li/MolecularMC/candida"
    #getGrAndConfig(chrRng,target_file_path,50000,25000)
    
    
    
    ### copy old result if needed
    #for Chromosome in ['1','2','3','4','5','6','7','R']: 
    #    os.system("cp ../v101/chr%s/*.finres chr%s"%(Chromosome,Chromosome))
    #print("copy done!")
    
    
    
    ### submit tasks for all sections to ITP cluster
    ### for calculate the optimized potential
    ### change accuracy need to change here (namesubfix) and file "runISS.py"
    #namesubfix="finres01" 
    #submitAllMissing(namesubfix)
    
    ### or submit quest for one chromosome
    ### format: submitquest(Chromosome_name,range_of_section)
    #submitquest("R",[12])
    
    
    
    ### calculate conformations and compressibilities
    ### default box size range for compressibilities is np.arange(1,41,1), 
    ### it can be changed in function calCompressibility()
    ### default box size cutoff is 15 in function calCompressibility()
    #namesubfix="finres0.1" 
    
    #chrRng = ['R'] 
    #calConformation(namesubfix)
    #calCompressibility("compress01_100000.out")
    #os.system("rm */*.confor")
    
    
          