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
#from pylab import *

#####
"""
this is file to control and run all the function

change area length need to change this file (maxsetion, calCompress), 1getGr.cpp
"""
####

### set chromosome index
Chromosome = 'R'

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
if False:##need install noisyopt
    submitquest(Chromosome,range(maxsection))



if False:
    chrRng = ['1','2','3','4','5','6','7','R'] 
    for ichr in chrRng:
        os.chdir("chr%s"%ichr)
        isecrng = [int(f[5:-6]) for f in os.listdir() if f.endswith(".midgr")]
        #print(np.amax(isecrng))
        misarr = []
        for i in range(np.amax(isecrng)+1):
            try:
                data = np.genfromtxt("chr%s-%d.finres"%(ichr,i))
                #print(i,data)
            except:
                #print(i," is missing.")
                misarr.append(i)
        print(ichr,misarr)
        os.chdir("..")
        
        submitquest(ichr,misarr)

if False:
    os.system("g++ 5getConfor.cpp")
    #chrRng = ['7'] 
    chrRng = ['1','2','3','4','5','6','7','R'] 
    for ichr in chrRng:
        os.chdir("chr%s"%ichr)
        isecrng = [int(f[5:-7]) for f in os.listdir() if f.endswith(".finres")]
        isecrng.sort()
        #print(isecrng)
        for i in isecrng:
            os.system("../a.out %d %s"%(i,ichr))    
        os.chdir("..")

def calCompress(ichr,ita):
    dataname = "chr%s-%d.confor"%(ichr,ita)
    global boxrng

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
    #plot(75000/boxrng,stdarr,'*-',color=col[i])
    x = boxrng#/75000.
    y = np.array(rhoktKTarr)
    y = y[x>15]
    x = x[x>15]        
    A = np.vstack([x, np.ones(len(x))]).T
    sol = np.linalg.lstsq(A, y,rcond=None)
    v, c = sol[0]
    #print(i,v,c)
    return [i,v,c]+rhoktKTarr
if False:
    boxrng = np.arange(1,41,1)
    
    chrRng = ['1','2','3','4','5','6','7','R'] 
    #chrRng = ['2'] 
    for ichr in chrRng:
        os.chdir("chr%s"%ichr)
        isecrng = [int(f[5:-7]) for f in os.listdir() if f.endswith(".finres")]
        isecrng.sort()
        #print(isecrng)
        
        ffin = open("compress.out","w")
        ffin.write("0 0 0 ")
        for item in boxrng:
            ffin.write("%d "%item)
        ffin.write("\n")
            
        for i in isecrng:
            #print(ichr,i)
            irescomp = calCompress(ichr,i)
            print(ichr,irescomp[:3])
            for itemires in irescomp:
                ffin.write("%f "%itemires)
            ffin.write("\n")
        ffin.close() 
        os.chdir("..")
        
    
    
### plot compressibilities out as chromosome plot    
if False:    
    """
    Rough script to plot chromosome ideograms using data from UCSC
    from internet
    """
    
#    color_lookup = {
#                      '1': (1., 1., 1.),
#                    '2': (.6, .6, .6),
#                    '3': (.4, .4, .4),
#                    '4': (.2, .2, .2),
#                   '5': (0., 0., 0.),
#                      'acen': (.8, .4, .4),
#                      'gvar': (.8, .8, .8),
#                     'stalk': (.9, .9, .9),
#                   }
    def color_lookup(x):
        return plt.cm.Blues(x*8)
        #return plt.cm.autumn_r((np.clip(x,2,10)-2)/8.)
        
    height = 18
    spacing = 9
    
    def ideograms(fin):
        last_chrom = None
        xranges, colors = [], []
        ymin = 0
    
        for line in fin:
            chrom, start, stop, label, stain = line
            start = int(start)
            stop = int(stop)
            width = stop - start
            if chrom == last_chrom or (last_chrom is None):
                xranges.append((start, width))
                colors.append(color_lookup(stain))
                last_chrom = chrom
                continue
    
            ymin += height + spacing
            yrange = (ymin, height)
            yield xranges, yrange, colors, last_chrom
            xranges, colors = [], []
            xranges.append((start, width))
            colors.append(color_lookup(stain))
            last_chrom = chrom
    
        # last one
        ymin += height + spacing
        yrange = (ymin, height)
        yield xranges, yrange, colors, last_chrom
    
    fig = plt.figure()
    ax = fig.add_subplot(111)
    d = {}
    yticks = []
    yticklabels = []
    
    #chrRng = ['1','2','3','4','5','6','7'] 
    chrRng = ['R','7','6','5','4','3','2','1'] 
    mydata = []
    for ichr in chrRng:
        os.chdir("chr%s"%ichr)
        tcom = np.genfromtxt("compress.out",skip_header=1)
        for itemcom in tcom:
            i = itemcom[0]
            c = itemcom[2]
            mydata.append(["chr%s"%ichr,12500+i*50000,12500+(i+1)*50000,0,c])
        os.chdir("..")    
    #mydata = [["chr1",0,50000,0,1],
    #["chr1",50000,100000,0,2],
    #["chr1",100000,150000,0,1],["chr2",150000,200000,0,3],["chr2",200000,250000,0,4]]
    
    for xranges, yrange, colors, label in ideograms(mydata):
        coll = BrokenBarHCollection(xranges, yrange, facecolors=colors)
        ax.add_collection(coll)
        center = yrange[0] + yrange[1]/2.
        yticks.append(center)
        yticklabels.append(label)
        d[label] = xranges
    
    ax.axis('tight')
    ax.set_yticks(yticks)
    ax.set_yticklabels(yticklabels)
    ax.set_xticks([])
    plt.savefig("compress.pdf")
    plt.show()
    
### plot compressibilities out as curve        
if False:
    #isubplot=1
    
    chrRng = ['R','7','6','5','4','3','2','1'] 
    mydata = []
    for ichr in chrRng:
        os.chdir("chr%s"%ichr)
        tcom = np.genfromtxt("compress.out",skip_header=1)
        for itemcom in tcom:
            i = itemcom[0]
            c = itemcom[2]
            mydata.append(["chr%s"%ichr,12500+i*50000,12500+(i+1)*50000,0,c])
        mydata=np.asarray(mydata)[:,1:].astype(float)
        #plt.subplot(810+isubplot);isubplot+=1;
        #plt.figure()
        plt.plot((mydata[:,1]+mydata[:,2])/2,mydata[:,-1],"-*",alpha=0.6,label="chr%s"%ichr)
        #plt.title("chr%s"%ichr)
        mydata=[]
        os.chdir("..")     
    plt.legend()    
    plt.show()
    
    