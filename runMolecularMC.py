#! /remote/gpu05/anaconda3/bin/python3
# -*- coding: utf-8 -*-

#import sys
#sys.path.append('/home/li/bin/')from ctypes import cdll
from ctypes import *
#import numpy as np
#from pylab import *
#from scipy.optimize import curve_fit
from head import *
import pandas as pd
#### load the library lib = cdll.LoadLibrary('./cdll.so')
lib.c_loadTarget.argtypes = [c_char_p]lib.c_readConfig.argtypes = [c_char_p]lib.c_getBestGr.argtypes = [c_int]
lib.c_getBestGr.restype = c_double
lib.c_assignAndRun.argtypes = [c_double,c_double,c_double,c_double]
lib.c_assignAndRun.restype = c_double

#if False:
#    lib.c_SeedByTime()
#    xTarget,yTarget = np.genfromtxt("targetGr.gr").T
#    ngr = len(xTarget)
#    yBest = np.zeros(ngr)
#    yBestTmp = np.zeros(ngr)
#    dx = xTarget[1]-xTarget[0]
#
#    xlist = [2.5, 3.0, 3.1, 3.2, 3.3000000000000003, 3.4000000000000004, 3.5, 3.6, 3.7, 3.8000000000000003, 3.9000000000000004, 4.0, 4.1000000000000005, 4.2, 4.3, 4.4, 4.5, 4.6000000000000005, 4.7, 4.800000000000001, 4.9, 5.0, 5.1000000000000005, 5.2, 5.300000000000001, 5.4, 5.5, 5.6000000000000005, 5.7, 5.800000000000001, 5.9, 6.0, 6.1000000000000005, 6.2, 6.300000000000001, 6.4, 6.5, 6.6000000000000005, 6.7, 6.800000000000001, 6.9, 7.0, 7.1000000000000005, 7.2, 7.300000000000001, 7.4, 7.5, 7.6000000000000005, 7.7, 7.800000000000001, 7.9, 8.0, 8.1, 8.200000000000001, 8.3, 8.4, 8.5, 8.6, 8.700000000000001, 8.8, 8.9, 9.0, 9.1, 9.200000000000001, 9.3, 9.4, 9.5, 9.600000000000001, 9.700000000000001, 9.8, 9.9,30]
#    for x in xlist:
#        lib.c_xAdd(x)
#
#    oldp = 1e99
#    for i in range(10):
#        lib.c_yClear()
#        lib.c_yAdd(1e20)
#        for j in range(len(xlist)-2):
#            lib.c_yAdd(-10*np.random.rand()+5)
#        lib.c_yAdd(0)
#
#        if lib.c_xyEqual()==False:
#            print("size error")
#
#        lib.c_runPotInput(1,0.45,1e-6,0,800,200,500,1,500,ngr,dx)
#        for iBest in range(ngr):
#            yBestTmp[iBest] = lib.c_getBest(iBest)
#
#        newp = np.mean((yTarget-yBestTmp)**2)
#        if newp<=oldp:
#            yBest=yBestTmp
#            oldp=newp
#            print(newp)
#
#
#lib.c_SeedByTime()
#xTarget,yTarget = np.genfromtxt("targetGr.gr").T
#ngr = len(xTarget)
##yBest = np.zeros(ngr)
#yBestTmp = np.zeros(ngr)
#dx = xTarget[1]-xTarget[0]
#testcount=-1
#
#def interface(dum,*args):
#    global testcount,ngr,dx
#    testcount+=1
#    if testcount%100==0:
#        print(testcount,args[:10])
#
#    lib.c_yClear()
#    #lib.c_yAdd(1e20)
#    for j in range(len(args)):
#        lib.c_yAdd(args[j])
#    #lib.c_yAdd(0)
##    print(len(args))
##    lib.c_yPrint()
##    lib.c_xPrint()
##    print(lib.c_xyEqual())
##    exit()
#
#    lib.c_runPotInput(1,0.45,1e-6,0,800,200,500,1,500,ngr,dx)
#    for iBest in range(ngr):
#        yBestTmp[iBest] = lib.c_getBest(iBest)
#    return yBestTmp

#def initTarget(xlist):
#    curinp = interp1d(xTarget, -np.log(yTarget+1.1), kind='linear')
#    return curinp(xlist)
def polyN(xlist,*args):
    global testcount,ngr,dx
    testcount+=1
    if testcount%1==0:
        print(testcount,args[:10])

    lib.c_yClear()
    lib.c_yAdd(1e20)
    for j in range(1,len(xlist)-1):
        lib.c_yAdd(np.polyval(args,xlist[j]))
    lib.c_yAdd(0)
#    print(len(args))
#    lib.c_yPrint()
#    lib.c_xPrint()
#    print(lib.c_xyEqual())
#    exit()

    lib.c_runPotInput(1,0.45,1e-6,0,800,200,500,1,500,ngr,dx)
    for iBest in range(ngr):
        yBestTmp[iBest] = lib.c_getBest(iBest)
    return yBestTmp
def LJlike(r,epsilon,sigma,a,b):
    return 4*epsilon*((sigma/r)**a-(sigma/r)**b)
def LJlikeInterface(xlist,sigma):
    epsilon = 3
    a = 12
    b = 6
    global testcount,ngr,dx
    testcount+=1
    if testcount%1==0:
        print(testcount,epsilon,sigma,a,b)

    lib.c_yClear()
    #lib.c_yAdd(1e20)
    for j in range(len(xlist)):
        lib.c_yAdd(LJlike(xlist[j],epsilon,sigma,a,b))
    #lib.c_yAdd(0)
#    print(len(args))
#    lib.c_yPrint()
#    lib.c_xPrint()
#    print(lib.c_xyEqual())
#    exit()

    lib.c_runPotInput(1,0.45,1e-6,0,800,200,500,1,500,ngr,dx)
    for iBest in range(ngr):
        yBestTmp[iBest] = lib.c_getBest(iBest)
    return yBestTmp

#if False:
#    xlist = [2.5, 3.0, 3.1, 3.2, 3.3000000000000003, 3.4000000000000004, 3.5, 3.6, 3.7, 3.8000000000000003, 3.9000000000000004, 4.0, 4.1000000000000005, 4.2, 4.3, 4.4, 4.5, 4.6000000000000005, 4.7, 4.800000000000001, 4.9, 5.0, 5.1000000000000005, 5.2, 5.300000000000001, 5.4, 5.5, 5.6000000000000005, 5.7, 5.800000000000001, 5.9, 6.0, 6.1000000000000005, 6.2, 6.300000000000001, 6.4, 6.5, 6.6000000000000005, 6.7, 6.800000000000001, 6.9, 7.0, 7.1000000000000005, 7.2, 7.300000000000001, 7.4, 7.5, 7.6000000000000005, 7.7, 7.800000000000001, 7.9, 8.0, 8.1, 8.200000000000001, 8.3, 8.4, 8.5, 8.6, 8.700000000000001, 8.8, 8.9, 9.0, 9.1, 9.200000000000001, 9.3, 9.4, 9.5, 9.600000000000001, 9.700000000000001, 9.8, 9.9,30]
#    for x in xlist:
#        lib.c_xAdd(x)
#
#    #curve_fit(reCurve, xTarget, yTarget,p0=[20000,-20], maxfev=500000)
#
#    #params_0 = np.arange(1,11)
#    #params_0 = [0,0,1e1,1e2,1e3,1e3,1e4,1e4,1e4]
#    params_0 = [1]
#    print(len(params_0))
#
#    #popt, pcov = curve_fit(lambda x, *params_0: interface(x, *params_0),
#    #           xTarget, yTarget, p0=params_0,maxfev=100000)
#    popt, pcov = curve_fit(LJlikeInterface, xlist, yTarget, p0=params_0,maxfev=10000000)
#    print(popt)
#    #np.savetxt('savepotential1.pot', np.vstack((xlist,np.polyval(popt,xlist))).T,fmt='%e')

def polyInterface(xlist,*args):
    return np.polyval(args,xlist)
def fourInterface(xbin,*args):
    fity = args
    frng = np.linspace(-100,100,100)
    fftx = np.zeros(len(xbin))
    for t in range(len(xbin)):
        fftx[t] = np.sum(fity*exp(2*np.pi*np.complex(0,1)*xbin[t]*frng)).real 
    #fftx /= getSurface(xbin,fftx)
    return fftx
def fPeriodRestrict(x,begin,end):
    siz = end-begin
    if(x>end):
        x -= siz
    if(x<begin):
        x += siz
    return x
class classAnneal:    oldpara = []
    newpara = []
    vRndSize = []
    energy = 0
    digitFlag = True
    targetY = []
    smrat = 0    def __init__(self):
        pass
    def readyToRun(self,tag,cog,toldpara,tvRndSize):
        lib.c_loadTarget(str.encode(tag))
        lib.c_readConfig(str.encode(cog))
        lib.c_readyToRun()
        self.oldpara = copy.deepcopy(toldpara)
        self.newpara = copy.deepcopy(toldpara)
        self.energy = 99999
        self.vRndSize = tvRndSize
        self.smrat = 100## a smooth ratio of 100
        x,y = np.genfromtxt(tag).T
        self.targetY = smooth(y,self.smrat)
        
    def getCurEnergy(self):
        ycal = []
        for i in range(len(self.targetY)):
            ycal.append(lib.c_getBestGr(i))
        ycal2 = smooth(ycal,self.smrat)
        return np.mean((self.targetY-ycal2)**2)  
    def anneal(self,T,T_min,plambda,SAcount_max,fname):
        SAcount = 0
        i=0
        apFlag = True
        new_E = 99999
        ressave= []
           
        while(T>T_min and SAcount<SAcount_max):
            i = 0
            while(i<1):
                lib.c_assignAndRun(*self.oldpara)
                self.energy = self.getCurEnergy()
                self.updateNew()
                lib.c_assignAndRun(*self.newpara)
                new_E = self.getCurEnergy()
                apFlag = self.acceptOrNot(self.energy,new_E,T)
                if(apFlag):
                    self.oldpara[:] = self.newpara[:]
                    self.energy = new_E
                i+=1
                SAcount+=1
            T *= plambda
            print(SAcount,T,self.oldpara,self.energy)
            ressave.append([T,self.oldpara[0],self.oldpara[1],self.oldpara[2],
                self.oldpara[3],self.energy])
        np.savetxt(fname,ressave)
        return self.oldpara
    def acceptOrNot(self,oldp,newp,T):
        if(newp<=oldp):
            return True
        if(np.random.rand()<np.exp((oldp-newp)/T)):
            return True
        return False
    def updateNew(self):
        if(self.digitFlag):
            if(self.vRndSize[0]==0):
                self.newpara[0]=self.oldpara[0]
            else:
                self.newpara[0]=fPeriodRestrict(
                self.oldpara[0]+np.round(np.random.rand()*2*self.vRndSize[0]-
                self.vRndSize[0]),100,350)
            if(self.vRndSize[1]==0):
                self.newpara[1]=self.oldpara[1]
            else:
                self.newpara[1]=fPeriodRestrict(
                self.oldpara[1]+np.round(np.random.rand()*2*self.vRndSize[1]-
                self.vRndSize[1]),1,30)
            if(self.vRndSize[2]==0):
                self.newpara[2]=self.oldpara[2]
            else:
                self.newpara[2]=fPeriodRestrict(
                self.oldpara[2]+np.round(np.random.rand()*2*self.vRndSize[2]-
                self.vRndSize[2]),1,20)
            if(self.vRndSize[3]==0):
                self.newpara[3]=self.oldpara[3]
            else:
                self.newpara[3]=fPeriodRestrict(
                self.oldpara[3]+np.round(np.random.rand()*2*self.vRndSize[3]-
                self.vRndSize[3]),1,20)   
        else:
            #print("test")
            if(self.vRndSize[0]==0):
                self.newpara[0]=self.oldpara[0]
            else:
                self.newpara[0]=fPeriodRestrict(
                self.oldpara[0]+(np.random.rand()*2.*self.vRndSize[0]-
                self.vRndSize[0]),100,350)
            if(self.vRndSize[1]==0):
                self.newpara[1]=self.oldpara[1]
            else:
                self.newpara[1]=fPeriodRestrict(
                self.oldpara[1]+(np.random.rand()*2.*self.vRndSize[1]-
                self.vRndSize[1]),1,30)
            if(self.vRndSize[2]==0):
                self.newpara[2]=self.oldpara[2]
            else:
                self.newpara[2]=fPeriodRestrict(
                self.oldpara[2]+(np.random.rand()*2.*self.vRndSize[2]-
                self.vRndSize[2]),1,20)
            if(self.vRndSize[3]==0):
                self.newpara[3]=self.oldpara[3]
            else:
                self.newpara[3]=fPeriodRestrict(
                self.oldpara[3]+(np.random.rand()*2.*self.vRndSize[3]-
                self.vRndSize[3]),1,20)     

    

def SAMethod1(a,b,tagnam,cogname,savenam):
    print(tagnam,"simulated annealing started:\n",pd.read_csv(cogname,index_col=None))
    #a,b = 12,6
    #for a in range(2,19):
    #    for b in range(1,a):
    if True:    
        if True:
            ann = classAnneal()
            told = [180,10,a,b]
            tsiz = [30,3,0,0]
            ann.readyToRun(tagnam,cogname,told,tsiz)
            ann.digitFlag = False
            tmpPara = ann.anneal(100,0.001,0.993,1500,savenam)
                    
#if False:
#    ###seed
#    itarea="0"
#    grname="chr2-%s.LJgr"%itarea
#    lowname="chr2-%s.lowconfig"%itarea
#    midname="chr2-%s.midconfig"%itarea
#    highname="chr2-%s.highconfig"%itarea
#    
#    x,y = np.genfromtxt(grname).T
#    ind = np.where(x==500)[0][0]
#    x = x[:ind]
#    y = y[:ind]
#    y2 = smooth(y,100)
#    
#    lib.c_loadTarget(str.encode(grname))
#    lib.c_readConfig(str.encode(lowname))
#    lib.c_readyToRun()
#    
#    data = np.genfromtxt("chr2-0LJES10/pot.can").astype(float)
#    res = []
#    for idata in data:
#        print(idata)
#        energy = lib.c_assignAndRun(*idata[1:5])
#        exit()
#        #energy = lib.c_assignAndRun(163.0,2.0,15.,5.)
#        ycal = []
#        for i in range(len(x)):
#            ycal.append(lib.c_getBestGr(i))
#            ycal = ycal[:ind]
#        
#        ycal2 = smooth(ycal,100)
#        res.append(np.mean((y2-ycal2)**2))
#        print(res[-1])
#    res=np.asarray(res)
#    print(np.amin(res))
#    #plot(x,y2,'.-',alpha=0.3)
#    #plot(x,ycal2,'.-',alpha=0.3)
#    #show()
#    exit()
#
#    params_0=np.random.rand(200)
#    print(len(params_0))
#    popt, pcov = curve_fit(polyInterface, x, y, p0=params_0,maxfev=10000000)
#    print(popt)
#    plot(x,y,'.')
#    plot(x,polyInterface(x,*popt),'.')
#    ylim(0,8)
#    show()
def abIndexing(pos):
    abcount = 0
    for a in range(2,19):
        for b in range(1,a):
            #print(abcount,a,b)
            if abcount == abindex:
                return a,b
            abcount += 1
if True:
    lib.c_SeedByTime()
    #if sys.argv[1] == "dPair":
    #itarea="0"
    itarea = str(sys.argv[1])
    abindex = int(sys.argv[2])
    outer_define1,outer_define2 = abIndexing(abindex)  # 0-152
    
    print("a,b:",outer_define1,outer_define2)

    grname="chr2-%s.midgr"%itarea###midgr or LJgr or ...
    lowname="chr2-%s.lowconfig"%itarea
    midname="chr2-%s.midconfig"%itarea
    highname="chr2-%s.highconfig"%itarea
    savnam="pypotSAM1f%d_%d.midpot"%(int(outer_define1),int(outer_define2))
    SAMethod1(outer_define1,outer_define2,grname,midname,savnam)  
    
