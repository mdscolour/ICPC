#! /remote/gpu05/anaconda3/bin/python3
# -*- coding: utf-8 -*-
import os
import sys
import copy
import time
import numpy as np
#from head import *
#import pandas as pd
#from sklearn.cluster import KMeans
#import sklearn.cluster.k_means_

def bedSep(bedfilename="./ExpNum.bed"):
#if sys.argv[1] == "bedSep":
    global chrRng    

    i=0    
    for iPath in range(1):
        #chrRng = ['2L','2R','3L','3R','X','4']
        #chrRng = ['Y']
        #chrRng = list(range(1,20))+['X','Y'] 
        #chrRng = ['1','2','3','4','5','6','7','8','9','10',
        #'11','12','13','14','15','16','17'
        #,'18','19','20','21','22','X','Y']                          ## p3
        for chrIndex in chrRng:
            ######## Roman numeral or not        
            #nchr = 'chr'+int_to_roman(chrIndex) # may be capital
            nchr = 'chr'+str(chrIndex)
            
            wfile = open('_chr%s.like_bed'%(chrIndex),'w')
            wfile.write('chromosome start-site end-site size\n')
            wfile.write('Seperated for each chromosome\n')
            amx=0
            
            with open(bedfilename)as f:      ## name
                for line in f:
                    dataline = line.strip().split()
                    
                    #content.append(dataline)
                    if dataline[0] == nchr or dataline[0]==('Chr'+str(chrIndex)):
                        if int(dataline[2])>amx:
                            amx=int(dataline[2])
                        i += 1
                        wfile.write("%s\t%s\t%s\t%d\n"%(nchr,dataline[1],dataline[2],int(dataline[2])-int(dataline[1])))
                        #wfile.write("%s %s %d\n"%(dataline[1],dataline[2],int(dataline[3])))
                        #if int(dataline[3]) != int(dataline[2])-int(dataline[1]):
                        #    print "error"
                    #if i==10:
                    #    break          
            wfile.close()
    print("'bedSep' success. Read entry:",i," Maximum read position:",amx)
                
if __name__ == '__main__':
    ### global variable
    #chrRng = ['2']
    chrRng = ['1','2','3','4','5','6','7','R'] ###target chromosome index
    
    ###target file, when only one .bed file exists, excute:
    os.system("mv *.bed ExpNum.bed") 
    

    bedSep() ### got _chrX.like_bed from ExpNum.bed




    

    
    
