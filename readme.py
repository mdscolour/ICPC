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

#####
"""
this is file to control and run all the function
See "__main__" part
"""

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
### changing computational accuracy need to change "candidates" in file "runV10.py".
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
    chrRng = ['1','2','3','4','5','6','7','R'] 
    for ichr in chrRng:
        os.chdir("chr%s"%ichr)
        isecrng = [int(f[5:-6]) for f in os.listdir() if f.endswith(".midgr")]
        #isecrng = [int(f[5:-9]) for f in os.listdir() if f.endswith(".finres01")]
        isecrng.sort()
        #print(isecrng)
        for i in isecrng:
            os.system("mv chr%s-%d.%s chr%s-%d.finrestmp"%(ichr,i,namesubfix,ichr,i))
            os.system("../a.out %d %s"%(i,ichr))  
            os.system("mv chr%s-%d.finrestmp chr%s-%d.%s"%(ichr,i,ichr,i,namesubfix))  
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
    boxrng = np.arange(1,41,1)
    
    chrRng = ['1','2','3','4','5','6','7','R'] 
    #chrRng = ['R'] 
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
        

    

### the following is for clustering and its plot
def LJlike(r,sigma,epsilon,a,b):
    return 4*epsilon*((sigma/r)**a-(sigma/r)**b)
def chr2colorplot(col,axt):
    ### special for chr. 2
    res=[]
    for ipat in range(44):
        try:
            pat = "chr2/chr2-%d.finres"%ipat
            data = np.genfromtxt(pat)
            res.append(np.ravel(data)[:4])
        except:
            print("error reading of data: ",ipat)
            pass
    res = np.asarray(res).astype(float)
    #res = res[:44,:]
    if len(res) != 44:
        print("length error")
    #res = res[:44]
    #print(len(col))
    #print(*col,sep=",")
    #colori = ['b','r','r','g','c','g','m','g','y','y','m',
    #'c','c','r','r','m','r','r','r','g','y','m',
    #'c','r','m','c','c','r','y','g','m','y',
    #'g','r','m','m','g','m','g','r','g','g','c','k']
    
    #plt.figure()
    #for i in range(len(col)):
    #    plt.plot([i]*100,np.linspace(0.3,1,100),color=col[i])
    #    plt.plot([i]*100,np.linspace(0,0.3,100),color=colori[i])
    #plt.hlines(0.3,0,len(col))
    #savefig("chr2clustercolor.pdf")
    ##show()
    
    #axt.figure(figsize=(8,6))
    x = np.linspace(140,300,100000)
    for j in range(len(col)):
        axt.plot(x,LJlike(x,*res[j,:]),color=col[j],label="session "+str(j))
    
    axt.set_ylim(-3.0,1)
    #axt.vlines(167,-10,1)
    #axt.vlines(147,-10,1)
    #legend(loc="lower right")
    axt.set_xlabel(r"Distance $r$")
    axt.set_ylabel(r"Potential $P(r)$")
    #savefig("chr2potentialcolor.pdf")
    #plt.show()  
### start of the plot
### color_lookup and ideograms defined twice with different details
if False:
    res=[]
    chrRng = ['R','7','6','5','4','3','2','1'] 
    for ichr in chrRng:
        os.chdir("chr%s"%ichr)
        isecrng = [int(f[5:-6]) for f in os.listdir() if f.endswith(".midgr")]
        isecrng.sort()
        print(ichr,isecrng)
        
        for ipat in isecrng:
            try:
                pat = "chr%s-%d.finres"%(ichr,ipat)
    #            with open(pat, "r") as infile:
    #                for i in (infile.readlines() [-1:]):
    #                    res.append(i.split(" ")[1:-1])
    #                    break
                data = np.genfromtxt(pat)
                res.append(np.ravel(data)[:4])
            except:
                print("error reading of chr. %s data %d "%(ichr,ipat))
                pass
        os.chdir("..")
            
    res = np.asarray(res).astype(float)
    print("length",len(res))
    #print(res)
    res = res[:,:]
    #exit()
    
    x = np.linspace(140,300,10000)
    sample = []
    for j in range(len(res)):
        sample.append(LJlike(x,*res[j,:]))  
    sample = np.asarray(sample)#.reshape(-1, 1) 
          
    #cluster = AgglomerativeClustering(n_clusters=3, affinity='euclidean', linkage='ward')
    cluster = AgglomerativeClustering(n_clusters=4, affinity='euclidean', linkage='ward')
    clslvl = cluster.fit_predict(sample)
    #print(np.where(clslvl==0))
    #print(np.where(clslvl==1))
    #print(np.where(clslvl==2))
    #print(np.where(clslvl==3))
    #print(np.where(clslvl==4))
    #print(clslvl)
    #exit()


    color_lookup = {
                      '0': 'r',
                      '1': 'orange',
                    '2': 'g',
                    '3': 'c',
                    '4': 'm',
                   }
        
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
                xranges.append((start+5000, width-5000))
                colors.append(color_lookup[stain])
                last_chrom = chrom
                continue
    
            ymin += height + spacing
            yrange = (ymin, height)
            yield xranges, yrange, colors, last_chrom
            xranges, colors = [], []
            xranges.append((start+5000, width-5000))
            colors.append(color_lookup[stain])
            last_chrom = chrom
    
        # last one
        ymin += height + spacing
        yrange = (ymin, height)
        yield xranges, yrange, colors, last_chrom
    
    #fig = plt.figure()
    #ax = fig.add_subplot(111)
    figall, axsall = plt.subplots(2,2,figsize=(12,12))
    axsall = axsall.flatten()
    ax = axsall[1]
    axsall[0].text(-0.15, 1.05, "A", transform=axsall[0].transAxes, size=15, weight='bold')
    axsall[1].text(-0.15, 1.05, "B", transform=axsall[1].transAxes, size=15, weight='bold')
    axsall[2].text(-0.15, 1.05, "C", transform=axsall[2].transAxes, size=15, weight='bold')
    axsall[3].text(-0.15, 1.05, "D", transform=axsall[3].transAxes, size=15, weight='bold')
    d = {}
    yticks = []
    yticklabels = []
    
    
    mydata = []
    for ichr in chrRng:
        os.chdir("chr%s"%ichr)
        isecrng = [int(f[5:-6]) for f in os.listdir() if f.endswith(".midgr")]
        isecrng.sort()
        for i in isecrng:
            mydata.append(["chr%s"%ichr,12500+i*50000,12500+(i+1)*50000,0,999])
        os.chdir("..")    
    mydata = np.asarray(mydata)
    mydata[:,-1] = clslvl[:]
    
    mydatachr2col = []
    for imydatachr2 in mydata[np.where(mydata[:,0]=='chr2'),-1][0]:
        mydatachr2col.append(color_lookup[imydatachr2])
    ### plot potential for chr. 2
    chr2colorplot(mydatachr2col,axsall[0])
    #plt.show()
    #exit()

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
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["bottom"].set_visible(False)
    ax.spines["left"].set_visible(False)

### block density method to get compressibility    
#if True: 
    #figall, axsall = plt.subplots(1,2,figsize=(16,8))
    #axsall[0].text(-0.15, 1.0, "A", transform=axsall[0].transAxes, size=20, weight='bold')
    #axsall[1].text(-0.15, 1.0, "B", transform=axsall[1].transAxes, size=20, weight='bold')
    
    data = np.genfromtxt("chr2/compress.out")
    print(data.shape)

    vdata = data[1:,2]
    xplot = data[0,3:]
    ind = xplot>15
    x = xplot[ind]
    data = data[1:,3:]
    for i in range(len(data)):
        y = data[i,:]
        y = y[ind]
        
        A = np.vstack([x, np.ones(len(x))]).T
        sol = np.linalg.lstsq(A, y,rcond=None)
        v, c = sol[0]
        #print(c,vdata[i])
        #print(i,v,c)
        axsall[2].plot(xplot[ind],data[i,:][ind],'.-',color=plt.cm.Blues(8*c))
        axsall[2].plot(np.arange(0,41,1),v*np.arange(0,41,1)+c,color=plt.cm.Blues(8*c),alpha=0.3)
#    print(data[1:,0],data[1:,2])
#    for i in range(43):
#        plot(data[i+1,0],data[i+1,2],"*",color=col[i])
    axsall[2].scatter([0]*len(vdata),vdata,c=plt.cm.Blues(8*vdata),
        marker='^',label="Reduced isothermal compressibility")
    axsall[2].set_xlim(-3,41)
    axsall[2].legend()
    #plt.colorbar()
    cmap = plt.cm.Blues
    norm = mpl.colors.Normalize(vmin=np.amin(vdata)*8-0.1, vmax=np.amax(vdata)*8+0.1)
    cb = plt.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=cmap), 
       orientation='vertical',ax=axsall[2])
    cb.set_ticks(np.arange(0.06,0.13,0.02)*8)
    cb.set_ticklabels(["0.06","0.08","0.10","0.12"])
    axsall[2].set_xlabel(r"Number of block $M_b$")
    axsall[2].set_ylabel(r"Compressibility of block $\chi_T(L,L_0)$")
    #plt.savefig("blockdensity.pdf")
    #plt.show()   
        
### plot compressibilities out as chromosome plot    
#if True:    
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
    
    #fig = plt.figure()
    ax = axsall[3]
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
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["bottom"].set_visible(False)
    ax.spines["left"].set_visible(False)
    #plt.tight_layout()
    #plt.savefig("compresspotentialcolor.pdf")
    plt.show()

if False:  
    import scipy.cluster.hierarchy as shc
    from sklearn.cluster import AgglomerativeClustering
    from pylab import *
        
    figall, axsall = plt.subplots(1,2,figsize=(12,6))
    axsall[0].text(-0.15, 1.05, "A", transform=axsall[0].transAxes, size=15, weight='bold')
    axsall[1].text(-0.15, 1.05, "B", transform=axsall[1].transAxes, size=15, weight='bold')
    
    res=[]
    for ipat in range(44):
        try:
            pat = "chr2/chr2-%d.finres"%ipat
#            with open(pat, "r") as infile:
#                for i in (infile.readlines() [-1:]):
#                    res.append(i.split(" ")[1:-1])
#                    break
            data = np.genfromtxt(pat)
            res.append(np.ravel(data)[:4])
        except:
            print("error reading of data: ",ipat)
            pass
    res = np.asarray(res).astype(float)
    print("length",len(res))
    #print(res)
    #res = res[:44,:]
    if len(res) != 44:
        print("length error")
        
    x = np.linspace(140,300,10000)
    sample = []
    for j in range(44):
        sample.append(LJlike(x,*res[j,:]))  
    sample = np.asarray(sample)#.reshape(-1, 1) 

def chrindex_lookup(xi):
    xicount = 0
    chrRng = ['R','7','6','5','4','3','2','1'] 
    for ichr in chrRng:
        os.chdir("chr%s"%ichr)
        isecrng = [int(f[5:-6]) for f in os.listdir() if f.endswith(".midgr")]
        isecrng.sort()
        os.chdir("..")
        for ipat in isecrng:
            if xicount == xi:
                return [ichr,ipat]
            xicount += 1
                    
if False:
    import scipy.cluster.hierarchy as shc
    from sklearn.cluster import AgglomerativeClustering
    from pylab import *
    
    figall, axsall = plt.subplots(2,2,figsize=(12,12))
    axsall=axsall.flatten()
    axsall[0].text(-0.15, 1.05, "A", transform=axsall[0].transAxes, size=15, weight='bold')
    axsall[1].text(-0.15, 1.05, "B", transform=axsall[1].transAxes, size=15, weight='bold')
    axsall[2].text(-0.15, 1.05, "C", transform=axsall[2].transAxes, size=15, weight='bold')
    axsall[3].text(-0.15, 1.05, "D", transform=axsall[3].transAxes, size=15, weight='bold')
    
    res=[]
    chrRng = ['R','7','6','5','4','3','2','1'] 
    for ichr in chrRng:
        os.chdir("chr%s"%ichr)
        isecrng = [int(f[5:-6]) for f in os.listdir() if f.endswith(".midgr")]
        isecrng.sort()
        print(ichr,isecrng)
        
        for ipat in isecrng:
            try:
                pat = "chr%s-%d.finres"%(ichr,ipat)
    #            with open(pat, "r") as infile:
    #                for i in (infile.readlines() [-1:]):
    #                    res.append(i.split(" ")[1:-1])
    #                    break
                data = np.genfromtxt(pat)
                res.append(np.ravel(data)[:4])
            except:
                print("error reading of chr. %s data %d "%(ichr,ipat))
                pass
        os.chdir("..")
            
    res = np.asarray(res).astype(float)
    print("length",len(res))
    #print(res)
    res = res[:,:]
    #exit()
    
    x = np.linspace(140,300,10000)
    sample = []
    for j in range(len(res)):
        sample.append(LJlike(x,*res[j,:]))  
    sample = np.asarray(sample)#.reshape(-1, 1) 
if False:
    dend = shc.dendrogram(shc.linkage(sample, method='ward',metric='euclidean'),orientation="left",
      ax=axsall[0],color_threshold=500,no_labels=True)
    #print(dend['leaves'])
    #print(dend['color_list'])
    #print(dend['leaves_color_list'])
    #plt.show()
    #exit()
    
    indmaxarg = np.asarray(dend['ivl']).astype(int)
    #indmaxarg = np.flip(indmaxarg)
    mat = np.zeros((len(res),len(res)))
    #x = np.linspace(140,300,10000)
    for i in range(len(res)):
        for j in range(len(res)):
            dis = LJlike(x,*res[indmaxarg[i],:]) - LJlike(x,*res[indmaxarg[j],:])
            #print((dis**2).shape)
            #exit()
            mat[i,j] = np.sqrt(np.sum(dis**2))
            #mat[i,j] = np.mean(dis**2)
    axsall[1].pcolor(mat,cmap="coolwarm")
    #plt.colorbar(ax=axsall[1])
    cmap = plt.cm.coolwarm
    norm = mpl.colors.Normalize(vmin=np.amin(mat), vmax=np.amax(mat))
    cb = plt.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=cmap), 
       orientation='vertical',ax=axsall[1])

    summat = np.sum(mat,axis=0)
    for isample in sample:
        axsall[2].plot(x,isample,'c')
    axsall[2].plot(x,sample[indmaxarg[np.argmax(summat)]],'m',label="lowest similarity")
    axsall[2].plot(x,sample[indmaxarg[np.argmin(summat)]],'r',label="highest similarity")
    axsall[2].set_xlabel("")
    axsall[2].legend()
    
    
#    print(chrindex_lookup(indmaxarg[np.argmax(summat)]))
#    print(chrindex_lookup(indmaxarg[np.argpartition(summat,-2)[-2]]))
#    print(chrindex_lookup(indmaxarg[np.argmin(summat)]))
#    print(np.argmax(summat))
#    print(np.argpartition(summat,-2)[-2])
#    print(np.argmin(summat))
    
    print(chrindex_lookup(275))
    print(chrindex_lookup(276))
    print(chrindex_lookup(277))
    print(chrindex_lookup(278))
    print(chrindex_lookup(234))
    print((np.where(indmaxarg==275)))
    print((np.where(indmaxarg==276)))
    print((np.where(indmaxarg==277)))
    print((np.where(indmaxarg==278)))
    print((np.where(indmaxarg==234))) 
    #exit()
    
#    for isample in sample:
#        axsall[3].plot(x,isample,'c')
#    axsall[3].plot(x,sample[indmaxarg[2]],'b',label="class 1")
#    axsall[3].plot(x,sample[indmaxarg[50]],'k',label="class 2")
#    axsall[3].plot(x,sample[indmaxarg[150]],'peru',label="class 3")
#    axsall[3].legend()
    axsall[3].axis("off")
    
    #plt.tight_layout()
    plt.savefig("dendro2color.pdf")
    plt.show()     
if False:
    fig, ax = plt.subplots(figsize=(12,12))
    data = np.genfromtxt("chr2/compress.out")
    vdata = data[1:,2].reshape(1,-1)
    
    plt.pcolor(vdata*8,cmap=plt.cm.Blues)
    #ax.axis("off")
    ax.get_yaxis().set_visible(False)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["bottom"].set_visible(False)
    ax.spines["left"].set_visible(False)
    plt.xticks(np.arange(0,44,5),np.arange(0,44,5)*50000+12500)
    #plt.tight_layout()
    plt.savefig("chr2compcolor.svg")
    #plt.show()
if False:
    target_file = "/scratch/li/MolecularMC/candida/_chr%s.like_bed"%(2)
    data = np.genfromtxt(target_file,skip_header=2,usecols=(1,2)).astype(int)

    st = 12500+50000*8+5600
    ed = st+4000
    lplot = ed-st
    
    res = np.zeros(lplot)
    for ust,ued in data:
        if ued<st:
            continue
        if ust>ed:
            print(st,ed,ust)
            break
        res[max(ust-st,0):min(ued-st,lplot-1)] += 1
    plt.figure(figsize=(12,3))
    plt.plot(range(st,ed),res,"r")
    #plt.ylabel("Nucleosome position")
    #plt.xlabel("Position (bp)")
    plt.savefig("chr2inps.svg")
    plt.show()
    
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
    ### change accuracy need to change here (namesubfix) and file "runV10.py"
    #namesubfix="finres01" 
    #submitAllMissing(namesubfix)
    
    ### or submit quest for one chromosome
    ### format: submitquest(Chromosome_name,range_of_section)
    #submitquest("R",[12])
    
    
    
    ### calculate conformations and compressibilities
    ### default box size range for compressibilities is np.arange(1,41,1), 
    ### it can be changed in function calCompressibility()
    ### default box size cutoff is 15 in function calCompressibility()
    namesubfix="finres01" 
    #calConformation(namesubfix)
    calCompressibility("compress%s.out"%(namesubfix[-2:]))
    os.system("rm */*.confor")
    
          