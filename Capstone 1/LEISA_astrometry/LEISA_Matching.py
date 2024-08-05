#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov  2 11:32:23 2020

@author: Dennis Houlihan

Explanation: This code is used to match objects from a star catalog with the 
             calculated J2000 RA and Dec coordinate of each pixel in an LEISA
             image. This is specific to LEISA data as RA_DEC function uses extensions 
             that are specific to LEISA FITS images.
"""
from LEISA_Astrometry import RA_DEC
import numpy as np
from numpy import zeros
import pandas as pd
import glob
from astropy.visualization import ZScaleInterval 
from astropy.io import fits
import matplotlib.pyplot as plt
import scipy.optimize as opt;
np.warnings.filterwarnings('ignore')

#%%
files =sorted(glob.glob(('/Users/dennishoulihan/Documents/Capstone/Plutocruise-data/*.fit')))
Gaia = pd.read_csv(r'/Users/dennishoulihan/Documents/Capstone/Plutocruise-data/lsb_304c_Gaiabright.csv')
Gaia = np.array(pd.DataFrame(Gaia))

image = files[0]

def Match(image,catalog):
    """ Match indices of catalog objects to pixels in an image based on distance.
    
    Other Functions Used
    --------------------
    RA_DEC : Function in LEISA_Astrometry.py which calculates RA and Dec of each pixel
             in a LEISA image
             
    Parameters
    ----------
    image: a LEISA image 
    catalog: any star catalog
    
    Returns
    -------
    minimum: array shape (N,9) where N is number of objects in the catalog
    
             first column has indices of objects in the catalog
             second column has row indices (x) of pixel with minimum distance from object
             third column has column indices (y) of pixel with minimum distance from object
             fourth column has the calculated distances
             fifth column has RA coordinate of object
             sixth column has dec coordinate of object
             seventh column has objects R mean mag
             eighth column has RA coordinate of matched pixel
             ninth column has dec coordinate of matched pixel
    """
    
    coordinates = RA_DEC(image)
    
    ra = coordinates[0][:]  # Right ascension Coordinates for the image
    dec = coordinates[1][:] # Declination Coordinates for the image
    
    distance = zeros((len(catalog),256,256))    
    minimum = zeros((len(catalog),9))               # Creating empty arrays 
    index_min = zeros((len(catalog),3),dtype = int)

    
    for k in range(len(catalog)):
        for i in range(len(ra[2])):
            for j in range(len(ra[2])):
    
                    # Calculating distance from each object (Distance formula)
                    distance[k,i,j] = (((ra[3,i,j] - catalog[k,1])*np.cos((np.pi/180)*dec[3,i,j]))**2 + (dec[3,i,j]-catalog[k,3])**2)**0.5
                    
        # Taking the x and y indices of the pixels closest to the objects
        index_min[k,:] = np.unravel_index(np.argmin(distance[k,:,:], axis=None), distance.shape)

        minimum[k,0] = k # corresponds to the gaia object indices
        minimum[k,1] = index_min[k,1] # corresponds to row indices of pixels
        minimum[k,2] = index_min[k,2] # corresponds to column indices of pixels
        minimum[k,3] = distance[k,index_min[k,1],index_min[k,2]] # Calculated distance
        minimum[k,4] = catalog[k,1]  # catalog object RA   **(make sure the indices are correct)**
        minimum[k,5] = catalog[k,3]  # catalog object Dec
        minimum[k,6] = catalog[k,8]  # catalog object R mean mag
        minimum[k,7] = ra[3,index_min[k,1],index_min[k,2]]  # Matched RA coordinates
        minimum[k,8] = dec[3,index_min[k,1],index_min[k,2]]   # Matched Dec coordinates

    return minimum


matched115 = Match(image, Gaia)
#%%
openf = fits.open(image)

data = openf[0].data[3]
#data = data[0]+data[1]+data[2]+data[3]
#data = data / openf[3].data
interval = ZScaleInterval() 
interval = interval.get_limits(data) # getting min/max values 
#%% Plotting image
plt.figure(figsize = (7,7))
plt.axis('off')
plt.imshow(data, cmap='gray',norm = None, interpolation='nearest', vmin = interval[0], vmax = interval[1])
#plt.scatter(matched115[:,1],matched115[:,2], marker = 'o',facecolors='none', edgecolors='r')
#plt.gca().invert_yaxis()
plt.savefig("image.png", bbox_inches='tight')
#plt.colorbar()

#%%   Saving data as a FITS File

hdu = fits.PrimaryHDU(data)
hdul = fits.HDUList([hdu])
#hdul.writeto('summed304.fits')








