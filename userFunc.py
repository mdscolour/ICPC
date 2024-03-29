#! /remote/gpu05/anaconda3/bin/python3
# -*- coding: utf-8 -*-

import os
import sys
import copy
import time
import numpy as np
from scipy.interpolate import interp1d

### preparation of .like_bed file, .bed file is needed
def quickBedSep():
    os.system("./bedSep.py") 
    
def getGrAndConfig(chrRng,chrLength,target_file_pathT,seclen=50000,inclen=25000):
    #global chrRng
    #global chrLength
    for iChromosome in range(len(chrRng)):  
        Chromosome = chrRng[iChromosome] 
        ### set maximum length
        maxlen = chrLength[iChromosome]
        
        maxsection = int(maxlen/seclen);
        if maxsection*seclen+inclen > maxlen:
            maxsection-=1;
       
       
        ### preparation of .lowconfig and .midgr file        
        ### create .midgrpre and .preconfig. Run in the main folder.
        target_file = "%s/_chr%s.like_bed"%(target_file_pathT,Chromosome)
        ### compile
        os.system("g++ -o 1getGr.out 1getGr.cpp")
        ### create folders if not exist
        if not os.path.exists("chr%s"%Chromosome):
            os.system("mkdir chr%s"%Chromosome)
        else:
            print("Folder chr%s exists, write in that folder."%Chromosome)
            
        ### run with 6 paramters
        ### the default Gr is calculated to 1000 with stepsize 5
        ### i.e. in range(0,1000,5). Stepsize 5 is according to data quality.
        ### 1000 is decided by the scale of interest. The change of 
        ### these two paremeters need to go into 1getGr.cpp in function getCanGr() 
        os.system("./1getGr.out %s %d %s chr%s %d %d"%
        (Chromosome,maxsection,target_file,Chromosome,seclen,inclen))
        
        ### from .preconfig to .lowconfig. Run in the subfolder.
        os.system("g++ -o 2getConfig.out 2getConfig.cpp")
        os.chdir("chr%s"%Chromosome)
        for ita in range(maxsection):
            os.system("../2getConfig.out %d %s"%(ita,Chromosome))
        #os.chdir("..")
        
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
        print("For chr. %s, creation of RDF & config file success."%Chromosome)

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
def calConformation(chrRng,namesubfix="finres",mcs=100000):
    os.system("g++ -o 5getConfor.out 5getConfor.cpp")
    #chrRng = ['R'] 
    #chrRng = ['1','2','3','4','5','6','7','R'] 
    #global chrRng
    
    for ichr in chrRng:
        os.chdir("chr%s"%ichr)
        isecrng = [int(f.split(".")[0].split("-")[1]
                  ) for f in os.listdir() if f.endswith(namesubfix)]
        isecrng.sort()

        for i in isecrng:
            ### third parameter is the conformation number, times 10 is MSC
            #os.system("mv chr%s-%d.%s chr%s-%d.finrestmp"%(ichr,i,namesubfix,ichr,i))
            os.system("../5getConfor.out %d %s %d confor %s"%(i,ichr,mcs,namesubfix))  
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
def calCompressibility(chrRng,nameout="compress.out"):
    #boxrng = np.arange(1,41,1)
    boxrng = np.arange(16,41,2)
    
    #global chrRng
    #chrRng = ['1','2','3','4','5','6','7','R'] 
    #chrRng = ['1'] 
    for ichr in chrRng:
        os.chdir("chr%s"%ichr)
        isecrng = [int(f[5:-7]) for f in os.listdir() if f.endswith(".confor")]
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
            #print(ichr,irescomp[:3])
            for itemires in irescomp:
                ffin.write("%f "%itemires)
            ffin.write("\n")
        ffin.close() 
        os.chdir("..")
        
### combine the data into one file    
def summarizeAllValue(chrRng,namesubfix,namecompress,
    namerespot="Potential01_100000.like_bed"):  
    res=[]
    #global chrRng
    #chrRng = ['R']#,'7','6','5','4','3','2','1'] 
    #chrRng = ['R','7','6','5','4','3','2','1'][::-1] 
    for ichr in chrRng:
        os.chdir("chr%s"%ichr)
        isecrng = [int(f.split(".")[0].split("-")[1]
                  ) for f in os.listdir() if f.endswith(namesubfix)]
        isecrng.sort()

        #tcom = np.genfromtxt(namecompress,skip_header=1)
        tcom = np.loadtxt(namecompress,skiprows=1,ndmin=2)
        
        for ipat in isecrng:
            cLoc = np.where(tcom[:,0]==ipat)
            if len(cLoc)==0:
                c = "n/a"
            else:
                c = tcom[cLoc[0][0],2]     
            
            try:
                #pat = "chr%s-%d.finres01"%(ichr,ipat)
                pat = "chr%s-%d.%s"%(ichr,ipat,namesubfix)
                data = np.genfromtxt(pat)
                res.append(["chr%s"%ichr,12500+50000*ipat,75000+50000*ipat,ipat]+np.ravel(data)[:4].tolist()+[c])
            except:
                print("error reading of chr. %s data %d "%(ichr,ipat))
                pass
        os.chdir("..")
            
    #res = np.asarray(res).astype(float)
    #print("length",len(res))
    #print(res)
    #pd.DataFrame(res).to_csv(namerespot,sep="\t",index=False,header=None)
    fout = open(namerespot,'w') 
    fout.write("chromosome start_site end_site section_index "+
    "sigma epsilon delta nu compressibility\n")
    fout.write("File storing Effective potential and compressibilities\n")
    for etRes in res:
        for etEtRes in etRes:
            if type(etEtRes) == str:
                fout.write("%s\t"%(etEtRes)) 
            elif type(etEtRes) == float or type(etEtRes) == np.float64:
                fout.write("%f\t"%(etEtRes))
            elif type(etEtRes) == int:
                fout.write("%d\t"%(etEtRes))
            else:
                print("Value error in:",etRes) 
        fout.write("\n") 
    fout.close() 
    
       
#if __name__ == '__main__':
    ### get iNPS data, here is only a function for chromosome separation
    #quickBedSep() ### separate _chrX.like_bed from .bed file
    
    
    
    ### get radial distribution function and config file   
    ### last two parameter decide the section length
    ### change section length need to change here
    ### target chromosome index
    #chrRng = ['1','2','3','4','5','6','7','R'] 
    ### target chromosome length
    #chrLength = [3188341,2231883,1799298,1603259,1190845,1033292,949511,2286237]
    #target_file_path = "/scratch/li/MolecularMC/candida"
    #getGrAndConfig(chrRng,chrLength,target_file_path,50000,25000)
    
    
    
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
    #calConformation(chrRng,namesubfix)
    #calCompressibility(chrRng,"compress01_100000.out")
    #os.system("rm */*.confor")
    
    
          