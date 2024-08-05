#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Oct 11 23:36:23 2020

@author: Dennis Houlihan

Explanation: This script uses quaternions given in the headers of LEISA images 
             to rotate the cartesian pointing vectors of each pixel into 
             J2000 pointing vectors. Then RA and Dec coordinates are calculated 
             from those resulting vectors.
"""
from glob import glob
from astropy.io import fits
import numpy as np
import math as mt
np.warnings.filterwarnings('ignore')

def mult(q1, q2):
    """ Multiply two quaternions

    Parameters
    ----------
    q1 : 4 element sequence
    q2 : 4 element sequence

    Returns
    -------
    q12 : array shape (4,) 
    """
    w1, x1, y1, z1 = q1
    w2, x2, y2, z2 = q2
    w = w1*w2 - x1*x2 - y1*y2 - z1*z2
    x = w1*x2 + x1*w2 + y1*z2 - z1*y2
    y = w1*y2 + y1*w2 + z1*x2 - x1*z2
    z = w1*z2 + z1*w2 + x1*y2 - y1*x2
    return np.array([w, x, y, z])

def conjugate(q):
    """ Conjugate of quaternion

    Parameters
    ----------
    q : 4 element sequence
       w, i, j, k of quaternion

    Returns
    -------
    conjq : array shape (4,)
       w, i, j, k of conjugate of `q`
    """
    return np.array(q) * np.array([1.0, -1, -1, -1])

def rotate_vector(v, q):
    """ Apply transformation in quaternion `q` to vector `v`

    Parameters
    ----------
    v : 3 element sequence
       3 dimensional vector
    q : 4 element sequence
       w, i, j, k of quaternion

    Returns
    -------
    vdash : array shape (3,)
       `v` rotated by quaternion `q`
       
    """
    varr = np.zeros((4,))
    varr[1:] = v
    return mult(q, mult(varr, conjugate(q)))[1:]

def RA_DEC(fitsfile):
    """ Rotate pointing vectors for each pixel into the J2000 reference frame
        and then convert into RA and Dec Coordinates.
        
    Parameters
    ----------
    fitsfile : a FITS image, specifically from LEISA New Horizons datasets 
    - reason being that it calls for extensions specific to those FITS files
    
    
    Returns
    -------
    coordinates: array shape (2,N, 256,256) where N is number of frames in image
                 first index (0,N,256,256) is right ascension
                 second index (1,N,256,256) is declination
                 
    """
    openf = fits.open(fitsfile) #opening a fitsfile
    header = openf[0].header # Getting the header

# There are 3 256x256 arrays corresponding to the X,Y,Z coordinates for the pointing vectors for each pixel    
    pointing_vX = openf[2].data[0] # Calling the pointing vectors (x,y,z)
    pointing_vY = openf[2].data[1]
    pointing_vZ = openf[2].data[2]
    
    quaternion = np.zeros((len(openf[7].data),4))
    J2000 = np.zeros((len(openf[7].data),3,256,256))
    ra = np.zeros((len(openf[7].data),256,256))            # Creating empty arrays
    dec = np.zeros((len(openf[7].data),256,256))
    vector = np.zeros((3,256,256))
    
    for k in range(len(openf[0].data)):
        quaternion[k] = np.array([openf[7].data[k,1],openf[7].data[k,2],openf[7].data[k,3],header['SPCQZ']])
        #quaternion[k] = np.array([header['SPCQA'],header['SPCQX'],header['SPCQY'],header['SPCQZ']])
        for i in range(len(pointing_vX)):
            for j in range(len(pointing_vY)):
                
                # The cartesian pointing vectors for each pixel given in the extension of the fits file
                vector[:,i,j] = np.array([pointing_vX[i,j],pointing_vY[i,j],pointing_vZ[i,j]])
                # Rotate the cartesian pointing vectors by the quaternion to get the J2000 pointing vector
                J2000[k,:,i,j] = rotate_vector(vector[:,i,j], quaternion[k])
                
                # Now that we have the J2000 pointing vectors, we can convert to RA and Dec coordinates 
                #   dec = arcsin(z)
                #   RA = 360 + arctan(y/x) 
                dec[k,i,j] = mt.degrees(np.arcsin(J2000[k,2,i,j]))
                ra[k,i,j] = 360 - abs(mt.degrees(np.arctan(J2000[k,1,i,j]/J2000[k,0,i,j])))
            
    coordinates = np.array([ra,dec])       
    return coordinates


files =sorted(glob('C:/Users/Dennis Houlihan/Documents/Capstone/Capstone 1/Plutocruise-data/*.fit'))

test = RA_DEC(files[4])
 


