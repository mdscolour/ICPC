import os
import sys
import copy
import time
#sys.path.append('/home/li/bin/')
from ctypes import cdll
from ctypes import *
import numpy as np
#from pylab import *
#from scipy.optimize import curve_fit
#import pandas as pd

snsec = sys.argv[1]
keystr = sys.argv[2]
ressep = sys.argv[3]
if keystr == "testRun":
    MCStrategyRng = [2,5,10,50,100,500,1000,5000,10000]
elif keystr == "normalRun":
    MCStrategyRng = [100,200,500,1000,5000,10000,20000,50000,100000]  
else:
    raise Exception("Wrong Keyword. No specific MCS list.")
    
def smooth(y, box_pts):
    box = np.ones(box_pts)/box_pts
    y_smooth = np.convolve(y, box, mode='same')
    return y_smooth

if(not os.path.isfile("./cdll.so")):
    os.system("g++ -fPIC -shared -o cdll.so cdll.cpp")

#### load the library
lib = cdll.LoadLibrary('./cdll.so')
lib.c_loadTarget.argtypes = [c_char_p]
lib.c_readConfig.argtypes = [c_char_p]
lib.c_getBestGr.argtypes = [c_int]
lib.c_getBestGr.restype = c_double
lib.c_assignAndRun.argtypes = [c_double,c_double,c_double,c_double]
lib.c_assignAndRun.restype = c_double
lib.c_setNumConfigs.argtypes = [c_int]
lib.c_getNumConfigs.restype = c_int

### prepare data files
tag = "%s.midgr"%snsec      
cog = "%s.lowconfig"%snsec 

### generating all candidates
### accuracy can be changed here!
#ressep = 1

if ressep == "1":
    ressep = float(ressep)
    # final potential file name
    resfile = "%s.finres1"%(snsec) 
    candidates = []
    for c in np.arange(2,21,1):
        for d in np.arange(max(1,c-3),c,1):
            for a in np.arange(140,171,1):
                for b in np.arange(1,17,1):
                    candidates.append([a,b,c,d])
    candidates = np.asarray(candidates)
else:
    ressep = float(ressep)
    # final potential file name
    resfile = "%s.finres%.1f"%(snsec,ressep) 
    # need low accuracy data
    oldresfile = "%s.finres1"%snsec 
    
    oldres = np.ravel(np.genfromtxt(oldresfile))[:4]
    candidates = []
    for c in np.arange(max(2,oldres[2]-5),oldres[2]+5,ressep):
        for d in np.arange(max(1,c-3),c,ressep):
            for a in [oldres[0]]:
                for b in np.arange(1,17,ressep):
                    candidates.append([a,b,c,d])
    candidates = np.asarray(candidates)


smrat = 100 ## a smoothing of 100 points, that is 1 in Gr curve
targetX,targetY = np.genfromtxt(tag).T
targetY = smooth(targetY,smrat)
targetN = len(targetX)

lib.c_loadTarget(str.encode(tag))
lib.c_readConfig(str.encode(cog))
lib.c_readyToRun()

def getCurMoment():
    global targetY
    ycal = []
    for i in range(targetN):
        ycal.append(lib.c_getBestGr(i))
    ycal2 = smooth(ycal,smrat)
    ycaldiv = targetY-ycal2
    return np.mean(ycaldiv),np.std(ycaldiv,ddof=1)
def getCurEnergy():
    global targetY
    ycal = []
    for i in range(targetN):
        ycal.append(lib.c_getBestGr(i))
    ycal2 = smooth(ycal,smrat)
    return np.mean((targetY-ycal2)**2)

### noisyopt start    
#objcount=0
#def paraToMoment(x):
#    global objcount
#    objcount+=1
#    if objcount%100==0:
#        print(objcount)
#    a = x[0]
#    b = x[1]
#    c = x[2]
#    d = x[3]
#    lib.c_assignAndRun(a,b,c,d)
#    return getCurMoment()
def obj(x):
    #global objcount
    #objcount+=1
    #if objcount%100==0:
    #    print(objcount)
    a = x[0]
    b = x[1]
    c = x[2]
    d = x[3]
    lib.c_assignAndRun(a,b,c,d)
    return getCurEnergy()
    
### calculate for all candidates in list "candidates"
### this part is suitable for parallel computing
### parallel computing needed to be adjusted individually
for pnumcon in MCStrategyRng:
#for pnumcon in [20000,50000,100000]:
    canN = len(candidates)
    print("iteration:%i,N:%i"%(pnumcon,canN))
    lib.c_setNumConfigs(pnumcon)
    
    canM1 = np.zeros(canN)
    for ican in range(canN):
        canM1[ican]= obj(candidates[ican])
    
    canN2 = int(0.25*canN)  ### selection ratio is 0.25
    if canN2 == 0:
        canN2 = 1
    argcan2 = canM1.argsort()[:canN2]
    candidates = candidates[argcan2]
    ### out of loop when only one remains
    if canN2 == 1: 
        break

print(argcan2)
canM1 = canM1[argcan2]
print(canM1)
print(candidates)
### if already outside the loop and multiple ramain, choose the best one
np.savetxt(resfile,candidates[0].reshape(1,4))

    
    
    
    