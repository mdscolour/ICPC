#! /remote/gpu05/anaconda3/bin/python3
# -*- coding: utf-8 -*-

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
tag = "chr2-%s.midgr"%snsec      
cog = "chr2-%s.lowconfig"%snsec    
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
def paraToMoment(x):
    global objcount
    objcount+=1
    if objcount%100==0:
        print(objcount)
    a = x[0]
    b = x[1]
    c = x[2]
    d = x[3]
    lib.c_assignAndRun(a,b,c,d)
    return getCurMoment()
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
    
candidates = []
for c in np.arange(2,21,1):
    for d in np.arange(1,c,1):
        for a in np.arange(140,171,1):
            for b in np.arange(1,15,1):
                candidates.append([a,b,c,d])
candidates = np.asarray(candidates)

for pnumcon in [100,200,500,1000,5000,10000,20000,50000,100000]:
    canN = len(candidates)
    print("iteration:%i,N:%i"%(pnumcon,canN))
    lib.c_setNumConfigs(pnumcon)
    
    canM1 = np.zeros(canN)
    for ican in range(canN):
        canM1[ican]= obj(candidates[ican])
    
    canN2 = int(0.3*canN)
    if canN2 == 0:
        canN2 = 1
    argcan2 = canM1.argsort()[:canN2]
    candidates = candidates[argcan2]
#print(argcan2)
canM1 = canM1[argcan2]
print(canM1)
print(candidates)
    
    
    
    