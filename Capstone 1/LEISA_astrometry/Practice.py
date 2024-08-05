#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov  2 11:32:23 2020

@author: dennishoulihan
"""
from LEISA_Astrometry import RA_DEC
import numpy as np
from numpy import zeros
import matplotlib.pyplot as plt
import scipy.optimize as opt;
import pandas as pd
from glob import glob
from astropy.io import fits
import os
import math as mt
np.warnings.filterwarnings('ignore')

Gaia = pd.read_csv(r'/Users/dennishoulihan/Documents/Stellar_Astrophysics/Bp_Rp12.csv')

files =sorted(glob('/Users/dennishoulihan/Documents/Capstone/Plutocruise-data/*.fit'))
coordinates = RA_DEC(files[0])

ra = coordinates[0]
dec = coordinates[1]

store = zeros((len(Gaia),1))            # Creating empty arrays needed later
storeindex = zeros((len(Gaia),3)) 
minimum = zeros((len(ra),3))

for k in range((len(Gaia))):
    for i in range(len(ra[2])):
        for j in range(len(dec[2])):

            # Calculating distance from each object (Distance formula)
                store[k,0] = (((ra[len(ra)/2,i,j] - Gaia[k,1])*np.cos(dec[len(dec)/2,i,j]))**2 + (dec[len(dec)/2,i,j]-Gaia[k,2])**2)**0.5
                storeindex[i,1] = i #remembering index
                storeindex[j,2] = j
        
   # minimum[i,0] = min(store)   #storing minimum values from store
    index_min = np.argmin(store) #stores indices
    minimum[k,0] = k #corresponds to the gaia object indices
    minimum[k,1] = storeindex[index_min,1] #corresponds to pixel indices
    minimum[k,2] = storeindex[index_min,2]

# minimum contains minumim distances and the corresponding values. 
# minimum[:,0] contains the indices that corresonds to BRIGHT (your array)
# minimum[:,1] contains the indices that corresonds to USNO (the catalog)
    
minimum[:,1:2] = minimum[:,1:2] + 1  # corrects index numbers (offset due to python starting at zero)





'''
# In this step, we are matching the rest of the columns from USNO with the index from minimum

othercolumns=np.zeros((len(minimum),8)) # make sure the shape is correct

for i in range(256):
    for j in range(len(minimum)):
        store=np.array_equiv(USNO[i,0],minimum[j,1]) # This matches rest of USNO
        if(store==True):                             # with the index # from minimum
            othercolumns[j,1:] = USNO[i,1:]

'''


















