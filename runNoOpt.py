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
tag = "chr2-0.LJgr"      
cog = "chr2-0.lowconfig"  
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
    if objcount%100==0:
        print(objcount)
    a = 162
    b = 3
    c = x[0]
    d = x[1]
    lib.c_assignAndRun(a,b,c,d)
    return getCurEnergy()
    #return getCurMoment()

#start_time = time.time()
#print(obj([12,6]))
#print("--- %s seconds ---" % (time.time() - start_time))

#lib.c_setNumConfigs(100)

res = minimizeCompass(obj, x0=[12, 6], deltatol=1, paired=False)
print(objcount)
print(res)
