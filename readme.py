#! /remote/gpu05/anaconda3/bin/python3
# -*- coding: utf-8 -*-

import os
import sys
import copy
import time
import numpy as np
from scipy.interpolate import interp1d
from pylab import *

#####
"""
this is file to control and run all the function

change area length need to change this file, 1getGr.cpp
"""
####

### set chromosome index
Chromosome = '1'

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

#maxlen = 3188341
maxsection = int(maxlen/50000);
if maxsection*50000+25000 > maxlen:
    maxsection-=1;

###preparation of .like_bed file
if False:
    os.system("./0bedSep.py") 

###preparation of .lowconfig and .midgr file
if False:
    target_file = "/scratch/li/MolecularMC/candida/_chr%s.like_bed"%Chromosome
    os.system("g++ 1getGr.cpp")
    try:
        os.system("mkdir chr%s"%Chromosome)
    except:
        pass
    os.system("./a.out %s %d %s chr%s"%
    (Chromosome,maxlen,target_file,Chromosome))
#if False:
    os.system("g++ 2getConfig.cpp")
    os.chdir("chr%s"%Chromosome)
    for ita in range(maxsection):
        os.system("../a.out %d %s"%(ita,Chromosome))
    #os.chdir("..")
#if False:
    ###count as step 3
    #os.chdir("chr%s"%Chromosome)
    patrng = np.asarray([f for f in os.listdir() if f.endswith("midgrpre")])
    patrng.sort()

    res = []
    xres = np.arange(0,600,0.01)## first peak <200, so choose maximum length 600 here
    for pat in patrng:
        x,y = np.genfromtxt(pat).T
        resinp = interp1d(x, y, kind='linear')
        #resinp = interp1d(x, y, kind='quadratic')
        yres = resinp(xres)
        np.savetxt(pat[:-3], np.vstack((xres,yres)).T,fmt='%f')
#if False:
    os.system("rm *pre*")
    os.chdir("..")
 
###run and get the potential, the .finres file   
if False:##need install noisyopt
    os.chdir("chr%s"%Chromosome)
    for i in range(maxsection):
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


if False:
    os.system("g++ 5getConfor.cpp")
    os.chdir("chr%s"%"1")
    #for ita in range(maxsection):
    for ita in range(43):
        os.system("../a.out %d 1"%ita)    
    
if False:
    col=['c', 'r', 'g', 'm', 'm', 'm', 'c', 'c', 'c', 'y', 'm', 'm', 'c', 'c', 'c', 'c', 'r', 'r', 'm', 'm', 'g', 'r', 'm', 'm', 'r', 'c', 'c', 'c', 'r', 'r', 'r', 'r', 'g', 'g', 'm', 'm', 'c', 'm', 'g', 'm', 'c', 'c', 'c']
    
    boxrng = np.arange(1,41,1)
    ffin = open("compress.out","w")
    ffin.write("0 0 0 ")
    for item in boxrng:
        ffin.write("%d "%item)
    ffin.write("\n")
    
    for i in range(43):
    #if True:
    #    i=0
        print(i)
        dataname = "chr1-%d.confor"%i

        #stdarr = []
        rhoktKTarr = []
        for boxsize in boxrng:
            #binshist = np.arange(0,75000-1+boxsize,boxsize)
            binshist = np.linspace(0,75000,boxsize+1)
            #binshistmid = (binshist[:-1]+binshist[1:])/2.0
            density = [] 
            with open(dataname)as f:            
                for line in f:
                    datline = np.asarray(line.split()).astype(float)
                    nsam = len(datline)
                    ### counting
                    ht = np.histogram(datline,binshist)
                    density.append(ht[0])
            print(boxsize,np.array(density).shape,np.sum(density[0]))
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
        #plot(75000/boxrng,stdarr,'*-',color=col[i])
        x = boxrng#/75000.
        y = np.array(rhoktKTarr)
        y = y[x>15]
        x = x[x>15]        
        A = np.vstack([x, np.ones(len(x))]).T
        sol = np.linalg.lstsq(A, y,rcond=None)
        v, c = sol[0]
        print(i,v,c)
        for item in [i,v,c]+rhoktKTarr:
            ffin.write("%f "%item)
        ffin.write("\n")
        plot(boxrng,rhoktKTarr,'*-',color=col[i],label="y=%.6fx+%.6f"%(v,c))
        fitx = np.linspace(np.amin(boxrng)-1,np.amax(boxrng)+1,100)
        plot(fitx,v*fitx+c,color="k",alpha=0.2)        
    xlabel(r"number of block $M_b$")
    ylabel(r"reduced isothermal compressibility $\chi_T(L,L_0)$")
    legend()
    #xscale("log")
    show()  
    ffin.close()  
    
    
    
    
    
    
    
    