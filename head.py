import numpy as np
#import matplotlib.pyplot as plt
#from mpl_toolkits.mplot3d import Axes3D
import sys
import time
from random import *

#from pylab import *
#from matplotlib import *
from scipy.interpolate import *

from scipy import stats
from scipy.signal import find_peaks_cwt

from scipy.optimize import curve_fit
from scipy.optimize import minimize

import scipy.ndimage as ndimage
import scipy

import os
import glob

from decimal import *
import ctypes

from numpy import NaN, Inf

import copy
import re
import warnings
import shutil

plt.style.use('seaborn-pastel')
plt.rcParams['lines.linewidth'] = 2

# peaks = scipy.signal.argrelextrema(smooth_df, np.greater_equal, order=100)[0] ##p3
def argpeakdet(v,delta, x = None):
    """
    Converted from MATLAB script at http://billauer.co.il/peakdet.html

    Returns two arrays

    function [maxtab, mintab]=peakdet(v, delta, x)
    %PEAKDET Detect peaks in a vector
    %        [MAXTAB, MINTAB] = PEAKDET(V, DELTA) finds the local
    %        maxima and minima ("peaks") in the vector V.
    %        MAXTAB and MINTAB consists of two columns. Column 1
    %        contains indices in V, and column 2 the found values.
    %
    %        With [MAXTAB, MINTAB] = PEAKDET(V, DELTA, X) the indices
    %        in MAXTAB and MINTAB are replaced with the corresponding
    %        X-values.
    %
    %        A point is considered a maximum peak if it has the maximal
    %        value, and was preceded (to the left) by a value lower by
    %        DELTA.

    % Eli Billauer, 3.4.05 (Explicitly not copyrighted).
    % This function is released to the public domain; Any use is allowed.

    """
    maxtab = []
    mintab = []

    if x is None:
        x = arange(len(v))

    v = asarray(v)

    if len(v) != len(x):
        sys.exit('Input vectors v and x must have same length')

    if not isscalar(delta):
        sys.exit('Input argument delta must be a scalar')

    if delta <= 0:
        sys.exit('Input argument delta must be positive')

    mn, mx = Inf, -Inf
    mnpos, mxpos = NaN, NaN

    lookformax = True

    for i in arange(len(v)):
        this = v[i]
        if this > mx:
            mx = this
            mxpos = x[i]
        if this < mn:
            mn = this
            mnpos = x[i]

        if lookformax:
            if this < mx-delta:
                maxtab.append(mxpos)
                mn = this
                mnpos = x[i]
                lookformax = False
        else:
            if this > mn+delta:
                mintab.append(mnpos)
                mx = this
                mxpos = x[i]
                lookformax = True

    return np.asarray(maxtab,dtype=int)
def peakdet(v, delta, x = None): #base on height
    """
    Converted from MATLAB script at http://billauer.co.il/peakdet.html

    Returns two arrays

    function [maxtab, mintab]=peakdet(v, delta, x)
    %PEAKDET Detect peaks in a vector
    %        [MAXTAB, MINTAB] = PEAKDET(V, DELTA) finds the local
    %        maxima and minima ("peaks") in the vector V.
    %        MAXTAB and MINTAB consists of two columns. Column 1
    %        contains indices in V, and column 2 the found values.
    %
    %        With [MAXTAB, MINTAB] = PEAKDET(V, DELTA, X) the indices
    %        in MAXTAB and MINTAB are replaced with the corresponding
    %        X-values.
    %
    %        A point is considered a maximum peak if it has the maximal
    %        value, and was preceded (to the left) by a value lower by
    %        DELTA.

    % Eli Billauer, 3.4.05 (Explicitly not copyrighted).
    % This function is released to the public domain; Any use is allowed.

    """
    maxtab = []
    mintab = []

    if x is None:
        x = arange(len(v))

    v = asarray(v)

    if len(v) != len(x):
        sys.exit('Input vectors v and x must have same length')

    if not isscalar(delta):
        sys.exit('Input argument delta must be a scalar')

    if delta <= 0:
        sys.exit('Input argument delta must be positive')

    mn, mx = Inf, -Inf
    mnpos, mxpos = NaN, NaN

    lookformax = True

    for i in arange(len(v)):
        this = v[i]
        if this > mx:
            mx = this
            mxpos = x[i]
        if this < mn:
            mn = this
            mnpos = x[i]

        if lookformax:
            if this < mx-delta:
                maxtab.append((mxpos, mx))
                mn = this
                mnpos = x[i]
                lookformax = False
        else:
            if this > mn+delta:
                mintab.append((mnpos, mn))
                mx = this
                mxpos = x[i]
                lookformax = True

    return array(maxtab), array(mintab)

def slideWin(data,win):
#    win=1
    #pedata=np.ones(len(data)-win+1)
    pedata=copy.deepcopy(data)
    for i in range(len(data)-win+1):
        pedata[i]=np.mean(data[i:i+win])
    return pedata


def int_to_roman(input):
   """
   Convert an integer to Roman numerals.

   Examples:
   >>> int_to_roman(0)
   Traceback (most recent call last):
   ValueError: Argument must be between 1 and 3999

   >>> int_to_roman(-1)
   Traceback (most recent call last):
   ValueError: Argument must be between 1 and 3999

   >>> int_to_roman(1.5)
   Traceback (most recent call last):
   TypeError: expected integer, got <type 'float'>

   >>> for i in range(1, 21): print int_to_roman(i)
   ...
   I
   II
   III
   IV
   V
   VI
   VII
   VIII
   IX
   X
   XI
   XII
   XIII
   XIV
   XV
   XVI
   XVII
   XVIII
   XIX
   XX
   >>> print int_to_roman(2000)
   MM
   >>> print int_to_roman(1999)
   MCMXCIX
   """
   if type(input) != type(1):
      raise( TypeError, "expected integer, got %s" % type(input))
   if not 0 < input < 4000:
      raise( ValueError, "Argument must be between 1 and 3999")
   ints = (1000, 900,  500, 400, 100,  90, 50,  40, 10,  9,   5,  4,   1)
   nums = ('M',  'CM', 'D', 'CD','C', 'XC','L','XL','X','IX','V','IV','I')
   result = ""
   for i in range(len(ints)):
      count = int(input / ints[i])
      result += nums[i] * count
      input -= ints[i] * count
   return result



def roman_to_int(input):
   """
   Convert a roman numeral to an integer.

   >>> r = range(1, 4000)
   >>> nums = [int_to_roman(i) for i in r]
   >>> ints = [roman_to_int(n) for n in nums]
   >>> print r == ints
   1

   >>> roman_to_int('VVVIV')
   Traceback (most recent call last):
    ...
   ValueError: input is not a valid roman numeral: VVVIV
   >>> roman_to_int(1)
   Traceback (most recent call last):
    ...
   TypeError: expected string, got <type 'int'>
   >>> roman_to_int('a')
   Traceback (most recent call last):
    ...
   ValueError: input is not a valid roman numeral: A
   >>> roman_to_int('IL')
   Traceback (most recent call last):
    ...
   ValueError: input is not a valid roman numeral: IL
   """
   if type(input) != type(""):
      raise( TypeError, "expected string, got %s" % type(input))
   input = input.upper()
   nums = ['M', 'D', 'C', 'L', 'X', 'V', 'I']
   ints = [1000, 500, 100, 50,  10,  5,   1]
   places = []
   for c in input:
      if not c in nums:
         raise( ValueError, "input is not a valid roman numeral: %s" % input)
   for i in range(len(input)):
      c = input[i]
      value = ints[nums.index(c)]
      # If the next place holds a larger number, this value is negative.
      try:
         nextvalue = ints[nums.index(input[i +1])]
         if nextvalue > value:
            value *= -1
      except IndexError:
         # there is no next place.
         pass
      places.append(value)
   sum = 0
   for n in places: sum += n
   # Easiest test for validity...
   if int_to_roman(sum) == input:
      return sum
   else:
      raise( ValueError, 'input is not a valid roman numeral: %s' % input)
def write1DArray(data,nam,tp="%d \n",val='w'):
    wfile = open(nam,val)
    for dd in data:
        wfile.write(tp%(dd))
    wfile.close()

def write2DArray(data,nam,tp="%d %d \n",val='w'):
    wfile = open(nam,val)
    for dd in data:
        wfile.write(tp%(dd[0],dd[1]))
    wfile.close()

def leastSquare(x,y):
    A = np.vstack([x, np.ones(len(x))]).T
    sol = np.linalg.lstsq(A, y,rcond=None)
    v, c = sol[0]
    fitx = np.linspace(min(x),max(x),100)
    fity = v*fitx+c
    return v,c,sol[1],fitx,fity

#axbx = ax[pltInd].get_position()
#axtemp = fig.add_axes([axbx.x0+0.5*axbx.width, axbx.y0+0.07*axbx.height, axbx.width*0.45, axbx.height*0.33])
#axtemp.set_xticklabels([])
#axtemp.plot(x,[0]*len(x),"-")
#for tick in axtemp.xaxis.get_major_ticks():
#     tick.label.set_fontsize(7)


def getSurface(x,y):
    #######normalize
    #surface = 0
    #for i in range(1,len(xrng)):
    #    surface += (xrng[i]-xrng[i-1])*wei[i]
    #wei /= surface

    surface = 0
    for i in range(1,len(x)):
        surface += (x[i]-x[i-1])*y[i]
    return surface

def derivative(x,y):
##central difference
    #x = np.array([1,2,3], dtype=np.float)
    #y = np.asarray([1,2,6])
    x = np.asarray(x)
    y = np.asarray(y)

    z1 = np.hstack((y[0], y[:-1]))
    z2 = np.hstack((y[1:], y[-1]))

    dx1 = np.hstack((0, np.diff(x)))
    dx2 = np.hstack((np.diff(x), 0))

    return (z2-z1) / (dx2+dx1)

def smooth(y, box_pts):
	box = np.ones(box_pts)/box_pts
	y_smooth = np.convolve(y, box, mode='same')
	return y_smooth
def linearModel(x,a,b):
    return a*x+b
def stretchExp(x,beta,nterm,const):
    #print beta,nterm
    return nterm*(np.exp(-x**beta)+const)
def exponentialModel(x,a,b,c):
    return a*np.exp(b*x)+c
def powerlawModel(x,a,b,c):
    return a*x**b+c
def gaus(x,a,x0,sigma):
    #global apeak
    return a*exp(-(x-x0)**2/(2*sigma**2))
def generalModel(x,a,b,l,t):
    return a*(x**l)*np.exp(b*(x**t))



