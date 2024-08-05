# -*- coding: utf-8 -*-
"""
Created on Sun Mar 28 12:46:34 2021

@author: Dennis Houlihan

Explanation: This code is used to identify stars in a LEISA image, create a table of 
             star cutouts, and generate a PSF from those cutouts.
"""
from glob import glob
from astropy.io import fits
import matplotlib.pyplot as plt
from astropy.visualization import simple_norm
from Cap2_functions import stars, PSF_maker
from astropy.visualization import ZScaleInterval 


# Calling the FITS file and getting its data
files =sorted(glob('C:/Users/Dennis Houlihan/Documents/Capstone/Capstone 1/LEISA_Astrometry/Vega-test/*.fit'))
image = fits.open(files[0])
data = image[0].data


# creating a summed image with different frames where vega is spread out
data_new = data[40,:,:] + data[56,:,:] + data[71,:,:] + data[88,:,:] + data[103,:,:]

# these two lines are for plotting the image
interval = ZScaleInterval()   
z = interval.get_limits(data_new)  

# Identifying the stars and creating table of them
stars_tbl= stars(data_new,0.003)

stars_tbl.remove_rows([1,5])  # manually removing stars I saw were misidentified
print('New Table',stars_tbl)

# function that creates PSF
epsf, std, fitted_stars = PSF_maker(data_new, stars_tbl, 4, 12)
    
# Plotting the generated ePSF
plt.figure()
norm = simple_norm(epsf.data, 'log', percent=99.)
plt.imshow(epsf.data, norm=norm, origin='lower', cmap='viridis')
plt.colorbar()