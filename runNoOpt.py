#! /remote/gpu05/anaconda3/bin/python3
# -*- coding: utf-8 -*-

from noisyopt import minimizeCompass

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
tag = "%s.midgr"%snsec      
cog = "%s.lowconfig"%snsec  
resfile = "%s.finres"%snsec 
if len(sys.argv)>2:
    snx0 = [sys.argv[2],sys.argv[3],sys.argv[4],sys.argv[5]]
    snx0 = np.asarray(snx0).astype(float)
else:
    snx0 = np.genfromtxt(resfile)

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

### prepare globals
smrat = 100
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
objcount=0
def obj(x):
    global objcount
    objcount+=1
    #if objcount%1000==0:
    #    print(objcount)
    a = x[0]
    b = x[1]
    c = x[2]
    d = x[3]
    lib.c_assignAndRun(a,b,c,d)
    return getCurEnergy()
    #return getCurMoment()

#start_time = time.time()
#print(obj([12,6]))
#print("--- %s seconds ---" % (time.time() - start_time))

lib.c_setNumConfigs(100)

res = minimizeCompass(obj, x0=snx0, deltatol=1, paired=False)
print(objcount)
print(res.success)
print(snx0)
print(res.x)
np.savetxt(resfile, np.vstack((snx0,res.x)))


