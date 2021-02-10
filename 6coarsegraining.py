
# -*- coding: utf-8 -*-
import numpy as np
#from pylab import *
import copy
import os
#import pandas as pd
#from head import *

if False:
    patrng = [f for f in os.listdir('./') if f.endswith("pot")]
    print(len(patrng))
    arr = []
    for pat in patrng:
        p1,p2,p3,p4,E = np.genfromtxt(pat).T
        arr.append([p3,p4,E])
    arr = np.asarray(arr)
    print(arr.shape)
    print('minimum ',arr[np.argmin(arr[:,-1]),:])
    H, edges = np.histogramdd(arr, bins = (10,10))
    print(np.argmax(H))
    #print(H.shape, edges[0].size, edges[1].size, edges[2].size)
if True:
    patrng = [f for f in os.listdir('./') if f.endswith("pot")]
    pat = patrng[0]
    data = np.genfromtxt(pat)
    arr = data[:,:]

    arr = np.asarray(arr)
    print(arr.shape)
    print('minimum ',arr[np.argmin(arr[:,-1]),:])
    #H, edges = np.histogramdd(arr[:,:2], bins = (10,10))
    #print(np.argmax(H))
    #print(H.shape, edges[0].size, edges[1].size, edges[2].size)









