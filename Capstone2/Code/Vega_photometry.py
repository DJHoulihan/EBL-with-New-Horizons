# -*- coding: utf-8 -*-
"""
Created on Tue Feb  2 19:07:54 2021

@author: Dennis Houlihan
"""

import numpy as np
from photutils import datasets
from photutils import find_peaks
from photutils import detect_threshold
from glob import glob
from astropy.visualization import ZScaleInterval 
from astropy.io import fits
from astropy.table import Table
from astropy.table import Column
from astropy.stats import sigma_clipped_stats
from astropy.nddata import NDData
from photutils.psf import extract_stars
import matplotlib.pyplot as plt
from astropy.visualization import simple_norm
from photutils import EPSFBuilder, CircularAperture, aperture_photometry
import pylab as py
from scipy.interpolate import RegularGridInterpolator
import scipy.optimize as opt;



y = np.array([62,66,65,63,61,58,61,64,67,63,60,58,61,64,66,64])
x = np.array([15,19,20,24,31,37,44,50,57,61,64,66,69,16,19,22])


plt.scatter(x,y)

#%%     Finding how far up the y-direction Vega travels each frame


def line(x,m,b):
    return m*x + b
x = np.array([15,19,20,24,31,37,44,50,57,61,64,66,69])

t = np.linspace(0,25,13)
parameters, covariance = opt.curve_fit(line,t,x)
perr = np.sqrt(np.diag(covariance))

plt.plot(t,line(t,*parameters),color='orange',linewidth=2)
plt.scatter(t,x)
print(parameters[0],parameters[1])

plt.title("Vega's Slope up the X-direction")
plt.xlabel('Frame, Initial = 39')
plt.ylabel('X-position on Image')
slope = parameters[0]
#%%  Creating a window that will follow the star as it moves through the image
 
y = np.array([62,66,65,63,61,58,61,64,67,63,60,58,61,64,66,64])
x = np.array([15,19,20,24,31,37,44,50,57,61,64,66,69,16,19,22])

files =sorted(glob('C:/Users/Dennis Houlihan/Documents/Capstone/Capstone 1/LEISA_Astrometry/Vega-test/*.fit'))
image = fits.open(files[0])
data0 = image[0].data
data = data0[38:129,:,:]

middle = x[0]

wymin = min(y) - 8
wymax = max(y) + 8

wxmin = np.zeros(len(data),dtype = int)
wxmax = np.zeros(len(data), dtype = int)
wxmin[0] = int(middle-10)
wxmax[0] = int(middle+10)

window = np.zeros(data.shape)
window = window[:,wxmin[0]:wxmax[0],wymin:wymax]
window[0] = data[0,wxmin[0]:wxmax[0],wymin:wymax]

for i in range(len(data)-1):
   
    wxmin[i+1] = wxmin[i] + 1
    wxmax[i+1] = wxmax[i] + 1
    window[i+1] = data[i+1,wxmin[i+1]:wxmax[i+1], wymin:wymax]
     
    
i=35
interval = ZScaleInterval() 
z = interval.get_limits(window[i])
plt.imshow(window[i], cmap='gray',norm = None, vmin = z[0] , vmax = z[1])


nrows = 6
ncols = 6
fig, ax = plt.subplots(nrows=nrows, ncols=ncols, figsize=(20, 20), squeeze=True)
ax = ax.ravel()

for i in range(nrows*ncols):
    norm = simple_norm(window[i], 'log', percent=99.)
    ax[i].imshow(window[i],norm=norm, origin='lower', cmap='viridis')
          
#%%   Getting the psf weighted flux for each window
from PSF_photometry import peaks_table, star_cutouts, PSF

peaks = []
stars = []
for x in range(3):
    peaks.append(peaks_table(window[i,:,:],0.002))
    stars.append(star_cutouts(window[i,:,:],peaks[i],2,2))

#psf = PSF(stars[:])
    























